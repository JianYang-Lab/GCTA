#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import read_table, run_gcta
from plink import write_plink
from test_errors import expect_error

MARKERS = [
    [0, 1, 2, 0, 1, 2],
    [0, 0, 1, 1, 2, 2],
    [2, 2, 2, 0, 1, 0],
    [None, 1, None, 0, 2, 1],
]
SEXES = [1, 2, 1, 2, 0, 1]


def make_fixture(work_dir):
    prefix = Path(work_dir) / "input" / "x"
    write_plink(prefix, MARKERS, chromosomes=[23] * len(MARKERS))
    lines = prefix.with_suffix(".fam").read_text(encoding="utf-8").splitlines()
    prefix.with_suffix(".fam").write_text(
        "".join(
            f"{line.split()[0]}\t{line.split()[1]}\t0\t0\t{sex}\t-9\n"
            for line, sex in zip(lines, SEXES)
        ),
        encoding="utf-8",
    )
    return prefix


def expected_rows():
    expected = []
    for index, marker in enumerate(MARKERS):
        valid_n = 0
        male_valid_n = 0
        dosage = 0.0
        male_hom_alt_n = 0
        male_het_n = 0
        for value, sex in zip(marker, SEXES):
            if sex == 0 or value is None:
                continue
            valid_n += 1
            dosage += value
            if sex == 1:
                male_valid_n += 1
                if value == 2:
                    male_hom_alt_n += 1
                elif value == 1:
                    male_het_n += 1
        valid_alleles = 2 * valid_n - male_valid_n
        frequency = (dosage - male_hom_alt_n - 0.5 * male_het_n) / valid_alleles
        expected.append((index, frequency, valid_alleles))
    return expected


def check_output(gcta, work_dir, threads):
    prefix = make_fixture(work_dir)
    output = run_gcta(gcta, work_dir, "freq", ["--bfile", prefix, "--freqx"], threads=threads)
    rows = read_table(output.with_suffix(".frq"))
    assert rows[0] == ["CHR", "SNP", "POS", "A1", "A2", "AF", "NCHROBS"]
    assert len(rows) == len(MARKERS) + 1
    for fields, (index, frequency, valid_alleles) in zip(rows[1:], expected_rows()):
        assert fields[:5] == ["23", f"rs{index}", str(1000 + index), "A", "G"]
        assert math.isclose(float(fields[5]), frequency, rel_tol=1e-6, abs_tol=1e-6)
        assert fields[6] == str(valid_alleles)
    return output.with_suffix(".frq")


def main():
    parser = argparse.ArgumentParser(description="GCTA X-chromosome frequency regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("output", "no_x"))
    args = parser.parse_args()
    work = Path(args.work)

    if args.scenario == "output":
        single = check_output(args.gcta, work / "single", 1)
        multi = check_output(args.gcta, work / "multi", 4)
        assert single.read_bytes() == multi.read_bytes()
    else:
        prefix = Path(args.data).resolve() / "gcta_legacy" / "toy"
        expect_error(
            args.gcta,
            work,
            "no_x",
            ["--bfile", prefix, "--freqx", "--out", work / "result"],
        )
    print(f"Passed X-chromosome frequency scenario: {args.scenario}")


if __name__ == "__main__":
    main()
