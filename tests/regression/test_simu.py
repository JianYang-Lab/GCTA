#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink
from synthetic import family_markers
from test_errors import expect_error


def make_fixture(work_dir):
    markers, chromosomes = family_markers()
    prefix = Path(work_dir) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    qtl = Path(work_dir) / "qtl.txt"
    qtl.write_text("rs1 0.2\nrs2 -0.1\n", encoding="utf-8")
    return prefix, qtl


def check_qtl(gcta, work_dir, case_control, threads):
    prefix, qtl = make_fixture(work_dir)
    arguments = [
        "--bfile", prefix,
        "--simu-cc" if case_control else "--simu-qt",
    ]
    if case_control:
        arguments.extend((20, 20, "--simu-k", 0.5))
    else:
        arguments.extend(("--simu-rep", 2))
    arguments.extend(("--simu-hsq", 0.5, "--simu-causal-loci", qtl, "--simu-seed", 12345))
    output = run_gcta(gcta, work_dir, "simulated", arguments, threads=threads)
    par_rows = read_table(output.with_suffix(".par"))
    assert par_rows[0] == ["QTL", "RefAllele", "Frequency", "Effect"]
    assert [row[0] for row in par_rows[1:]] == ["rs1", "rs2"]
    assert [row[3] for row in par_rows[1:]] == ["0.2", "-0.1"]
    pheno_rows = read_table(output.with_suffix(".phen"))
    assert len(pheno_rows) == 128
    if case_control:
        values = [row[2] for row in pheno_rows]
        assert values.count("2") == 20
        assert values.count("1") == 20
        assert values.count("-9") == 88
    else:
        assert all(len(row) == 4 for row in pheno_rows)
        assert all(math.isfinite(float(value)) for row in pheno_rows for value in row[2:])
    return output.with_suffix(".phen"), output.with_suffix(".par")


def main():
    parser = argparse.ArgumentParser(description="GCTA simulation regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("qt", "cc", "bad_seed"))
    args = parser.parse_args()
    if args.scenario == "bad_seed":
        prefix, qtl = make_fixture(Path(args.work))
        expect_error(
            args.gcta,
            Path(args.work),
            "bad_seed",
            ["--bfile", prefix, "--simu-qt", "--simu-causal-loci", qtl,
             "--simu-seed", 100, "--out", Path(args.work) / "result"],
        )
    else:
        single = check_qtl(args.gcta, Path(args.work) / "single", args.scenario == "cc", 1)
        multi = check_qtl(args.gcta, Path(args.work) / "multi", args.scenario == "cc", 4)
        assert_tables_close(single[1], multi[1])
    print(f"Passed simulation scenario: {args.scenario}")


if __name__ == "__main__":
    main()
