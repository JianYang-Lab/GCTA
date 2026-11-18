#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import assert_grm_equal, read_table, run_gcta
from plink import write_plink
from test_grm_x import make_fixture as make_x_fixture
from test_ld_legacy import MARKERS


def assert_v2(output):
    log = output.with_suffix(".log").read_text(encoding="utf-8")
    assert "Accepted options:" not in log, "Command silently fell back to the legacy parser"


def assert_legacy(output):
    log = output.with_suffix(".log").read_text(encoding="utf-8")
    assert "Accepted options:" in log, "Command did not enter the legacy parser"


def check_pair(gcta, work, name, prefix, arguments, legacy_arguments=()):
    outputs = []
    for path, extra, threads in (("v2", arguments, 1), ("legacy", [*arguments, *legacy_arguments], 4)):
        output = run_gcta(
            gcta,
            work / path,
            name,
            ["--bfile", prefix, *extra],
            threads=threads,
        )
        (assert_v2 if path == "v2" else assert_legacy)(output)
        outputs.append(output)
    return outputs


def main():
    parser = argparse.ArgumentParser(description="GCTA dual-path dispatch and equivalence tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", default="dual_paths", choices=("dual_paths",))
    args = parser.parse_args()
    work = Path(args.work)
    prefix = work / "input" / "standard"
    write_plink(prefix, MARKERS)

    modern, legacy = check_pair(args.gcta, work, "grm", prefix, ["--make-grm"], ["--save-ram"])
    assert_grm_equal(modern, legacy)

    modern, legacy = check_pair(args.gcta, work, "freq", prefix, ["--freq"], ["--save-ram"])
    modern_rows = modern.with_suffix(".frq").read_text(encoding="utf-8").splitlines()[1:]
    legacy_rows = legacy.with_suffix(".freq").read_text(encoding="utf-8").splitlines()
    assert len(modern_rows) == len(legacy_rows) == len(MARKERS)
    for modern_row, legacy_row in zip(modern_rows, legacy_rows):
        modern_fields = modern_row.split()
        legacy_fields = legacy_row.split()
        assert modern_fields[1] == legacy_fields[0]
        assert modern_fields[3] == legacy_fields[1]
        assert float(modern_fields[5]) == float(legacy_fields[2])

    modern, legacy = check_pair(args.gcta, work, "bed", prefix, ["--make-bed"], ["--save-ram"])
    for suffix in (".bed", ".bim"):
        assert modern.with_suffix(suffix).read_bytes() == legacy.with_suffix(suffix).read_bytes()
    modern_fam = [line.split() for line in modern.with_suffix(".fam").read_text().splitlines()]
    legacy_fam = [line.split() for line in legacy.with_suffix(".fam").read_text().splitlines()]
    assert [row[:5] for row in modern_fam] == [row[:5] for row in legacy_fam]
    assert all(row[5] == "nan" for row in modern_fam)
    assert all(row[5] == "-9" for row in legacy_fam)

    x_prefix = make_x_fixture(work / "x", 23)
    x_fam = read_table(x_prefix.with_suffix(".fam"))
    x_prefix.with_suffix(".fam").write_text(
        "".join(
            f"{row[0]}\t{row[1]}\t0\t0\t{2 if row[4] == '0' else row[4]}\t-9\n"
            for row in x_fam
        ),
        encoding="utf-8",
    )
    x_fam = x_prefix.with_suffix(".fam").read_text(encoding="utf-8").splitlines()
    x_prefix.with_suffix(".fam").write_text(
        "".join(
            "\t".join([*line.split()[:4], "2" if line.split()[4] == "0" else line.split()[4], line.split()[5]]) + "\n"
            for line in x_fam
        ),
        encoding="utf-8",
    )
    modern, legacy = check_pair(args.gcta, work / "x", "grm", x_prefix, ["--make-grm-xchr"], ["--save-ram"])
    assert_grm_equal(modern, legacy)
    print("Passed v2/legacy dispatch assertions and polymorphic-input output comparisons")


if __name__ == "__main__":
    main()
