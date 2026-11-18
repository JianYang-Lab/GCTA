#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import shutil
import subprocess

SCENARIOS = (
    "missing_input",
    "bad_magic",
    "truncated_bed",
    "unwritable_output",
    "grm_short_write",
    "legacy_invalid_option",
    "ld_invalid_window",
    "ld_file_sparse",
)


def expect_error(gcta, work_dir, name, arguments):
    work_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    command = [str(Path(gcta).resolve()), *map(str, arguments), "--thread-num", "1"]
    try:
        result = subprocess.run(
            command,
            cwd=work_dir,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise AssertionError(f"GCTA did not fail within 30 seconds: {name}") from error
    (work_dir / f"{name}.console").write_text(result.stdout, encoding="utf-8")
    assert result.returncode != 0, f"Invalid input unexpectedly succeeded: {name}"
    assert "error" in result.stdout.lower(), f"Missing error diagnostic: {name}"


def prepare_bed(data_dir, work_dir, scenario):
    source = Path(data_dir) / "gcta_legacy" / "toy"
    prefix = Path(work_dir) / "input" / "toy"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".bim", ".fam"):
        shutil.copyfile(Path(f"{source}{suffix}"), Path(f"{prefix}{suffix}"))
    data = Path(f"{source}.bed").read_bytes()
    if scenario == "bad_magic":
        data = b"\x00" + data[1:]
    else:
        data = data[:4]
    Path(f"{prefix}.bed").write_bytes(data)
    return prefix


def main():
    parser = argparse.ArgumentParser(description="GCTA I/O error regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    args = parser.parse_args()
    work = Path(args.work).resolve()
    output = work / "result"
    if args.scenario == "missing_input":
        prefix = work / "absent" / "missing"
        arguments = ["--bfile", prefix, "--make-grm", "--out", output]
    elif args.scenario == "legacy_invalid_option":
        arguments = ["--reml", "--reml-alg", 3, "--out", output]
    elif args.scenario == "ld_invalid_window":
        prefix = Path(args.data) / "gcta_legacy" / "toy"
        arguments = ["--bfile", prefix, "--ld-matrix", "--ld-wind", "invalid", "--out", output]
    elif args.scenario == "ld_file_sparse":
        work.mkdir(parents=True, exist_ok=True)
        score = work / "sparse.score.ld"
        score.write_text(
            "SNP chr bp freq mean_rsq snp_num max_rsq ldscore\n"
            "rs0 1 1000 0.2 0.1 1 0.1 1.1\n"
            "rs1 1 20000 0.3 0.1 1 0.1 1.2\n"
            "rs2 1 40000 0.4 0.1 1 0.1 1.3\n",
            encoding="utf-8",
        )
        arguments = ["--ld-file", score, "--ld-score-region", 10, "--out", output]
    elif args.scenario == "unwritable_output":
        prefix = Path(args.data) / "gcta_legacy" / "toy"
        arguments = ["--bfile", prefix, "--make-grm", "--out", work / "absent" / "result"]
    elif args.scenario == "grm_short_write":
        if not Path("/dev/full").exists():
            print("SKIP: /dev/full is unavailable")
            return
        work.mkdir(parents=True, exist_ok=True)
        for suffix in (".grm.bin", ".grm.N.bin"):
            link = Path(f"{output}{suffix}")
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to("/dev/full")
        prefix = Path(args.data).resolve() / "gcta_legacy" / "test"
        arguments = ["--bfile", prefix, "--make-grm", "--out", output]
    else:
        prefix = prepare_bed(args.data, work, args.scenario)
        arguments = ["--bfile", prefix, "--make-grm", "--out", output]
    expect_error(args.gcta, work, args.scenario, arguments)
    if args.scenario == "legacy_invalid_option":
        console = (work / f"{args.scenario}.console").read_text(encoding="utf-8")
        assert "Accepted options:" in console, "The legacy option parser was not entered"
    elif args.scenario == "ld_file_sparse":
        console = (work / f"{args.scenario}.console").read_text(encoding="utf-8")
        assert "at least 3 SNPs" in console, "Missing sparse-region diagnostic"
    print(f"Passed I/O error regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
