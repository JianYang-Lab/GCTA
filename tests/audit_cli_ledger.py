#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


V2_PATTERN = re.compile(r"vector<string> supported_flagsV2\s*=\s*\{(.*?)\};", re.DOTALL)
LEGACY_PATTERN = re.compile(r'strcmp\(argv\[i\],\s*"(--[^"]+)"\)')
COVERED_EVIDENCE = {"independent_oracle", "format_decode", "error_path", "characterization"}


def parsed_flags(source):
    main = (source / "src/main.cpp").read_text(encoding="utf-8")
    match = V2_PATTERN.search(main)
    if not match:
        raise AssertionError("Cannot locate supported_flagsV2 in src/main.cpp")
    v2 = set(re.findall(r'"(--[^"]+)"', match.group(1)))
    legacy_source = (source / "legacy/option.cpp").read_text(encoding="utf-8")
    legacy = set(LEGACY_PATTERN.findall(legacy_source))
    return v2, legacy


def available_tests(build):
    result = subprocess.run(
        ["ctest", "--test-dir", str(build), "--show-only=json-v1"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    return {test["name"] for test in json.loads(result.stdout)["tests"]}


def audit(source, build):
    ledger_path = source / "tests/cli_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    v2, legacy = parsed_flags(source)
    expected = v2 | legacy
    entries = ledger["flags"]
    actual = set(entries)
    assert actual == expected, (
        f"Ledger/parser mismatch; missing={sorted(expected - actual)}, "
        f"stale={sorted(actual - expected)}"
    )
    assert len(v2) == ledger["m1"]["expected_v2_unique_flags"], "Unexpected v2 flag count"
    legacy_mains = set(ledger["m1"]["legacy_main_flags"])
    assert legacy_mains <= legacy, f"Unknown legacy main flags: {sorted(legacy_mains - legacy)}"

    errors = []
    tests = available_tests(build) if build else set()
    for flag in sorted(expected):
        entry = entries[flag]
        expected_parsers = sorted(
            parser for parser, flags in (("legacy", legacy), ("v2", v2)) if flag in flags
        )
        if sorted(entry["parsers"]) != expected_parsers:
            errors.append(f"{flag}: parser paths {entry['parsers']} != {expected_parsers}")
        m1 = flag in v2 or flag in legacy_mains
        expected_scope = "M1" if m1 else "M2"
        if entry["scope"] != expected_scope:
            errors.append(f"{flag}: scope {entry['scope']} != {expected_scope}")
        status = entry["status"]
        if m1 and status not in ("covered", "excluded"):
            errors.append(f"{flag}: M1 status is {status}")
        if not m1 and status not in ("covered", "excluded"):
            errors.append(f"{flag}: M2 is incomplete; status is {status}")
        if status == "covered":
            if not entry.get("tests"):
                errors.append(f"{flag}: covered without a ctest mapping")
            if entry.get("evidence") not in COVERED_EVIDENCE:
                errors.append(f"{flag}: invalid evidence {entry.get('evidence')}")
            if not entry.get("details"):
                errors.append(f"{flag}: missing observable-behavior details")
            if not entry.get("dispatch"):
                errors.append(f"{flag}: missing dispatch evidence")
            if not entry.get("threads"):
                errors.append(f"{flag}: missing threading evidence")
        if status == "excluded" and not entry.get("details"):
            errors.append(f"{flag}: exclusion lacks attempted paths and required conditions")
        if build:
            for test in entry.get("tests", []):
                if test not in tests:
                    errors.append(f"{flag}: mapped ctest does not exist: {test}")
    if errors:
        raise AssertionError("CLI ledger audit failed:\n" + "\n".join(errors))
    return ledger


def run_covered(build, ledger):
    mapped = sorted({
        test
        for entry in ledger["flags"].values()
        if entry["status"] == "covered"
        for test in entry["tests"]
    })
    assert mapped, "No covered tests mapped"
    expression = "^(" + "|".join(re.escape(test) for test in mapped) + ")$"
    result = subprocess.run(
        ["ctest", "--test-dir", str(build), "--output-on-failure", "-j4", "-R", expression],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    print(result.stdout, end="")
    if result.returncode:
        raise AssertionError(f"Covered ctest command failed with exit code {result.returncode}")
    if "Not Run" in result.stdout or "***Skipped" in result.stdout:
        raise AssertionError("A covered ctest was skipped or not run")


def main():
    parser = argparse.ArgumentParser(description="Audit and execute the M1 CLI coverage ledger")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    ledger = audit(args.source.resolve(), args.build.resolve())
    if args.run:
        run_covered(args.build.resolve(), ledger)
    counts = {status: 0 for status in ("covered", "deferred", "excluded")}
    for entry in ledger["flags"].values():
        counts[entry["status"]] += 1
    print(f"CLI ledger passed: {len(ledger['flags'])} flags; {counts}")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, KeyError, ValueError, subprocess.CalledProcessError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
