#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink

MARKERS = [
    [0, 0, 0, 0, 1, 1, 2, 2],
    [0, 0, 0, 0, 1, 1, 2, 2],
    [2, 2, 2, 2, 1, 1, 0, 0],
    [0, 1, 0, 1, 0, 1, 0, 1],
]


def r_squared(first, second):
    first_mean = sum(first) / len(first)
    second_mean = sum(second) / len(second)
    covariance = sum((a - first_mean) * (b - second_mean) for a, b in zip(first, second)) / (len(first) - 1)
    first_variance = sum((value - first_mean) ** 2 for value in first) / (len(first) - 1)
    second_variance = sum((value - second_mean) ** 2 for value in second) / (len(second) - 1)
    return covariance * covariance / (first_variance * second_variance)


def expected_rsq():
    return [
        [0.0 if first_index == second_index else r_squared(first, second)
         for second_index, second in enumerate(MARKERS)]
        for first_index, first in enumerate(MARKERS)
    ]


def check_score(gcta, work_dir, threads):
    prefix = Path(work_dir) / "input" / "ld"
    write_plink(prefix, MARKERS)
    output = run_gcta(
        gcta,
        work_dir,
        "score",
        ["--bfile", prefix, "--ld-score", "--ld-wind", 1000],
        threads=threads,
    )
    rows = [row for row in read_table(output.with_suffix(".score.ld")) if row]
    expected = expected_rsq()
    assert rows[0] == ["SNP", "chr", "bp", "MAF", "mean_rsq", "snp_num", "max_rsq", "ldscore"]
    for index, fields in enumerate(rows[1:]):
        observed = [value for value in MARKERS[index] if value is not None]
        frequency = sum(observed) / (2.0 * len(observed))
        pair_values = [expected[index][other] for other in range(len(MARKERS)) if other != index]
        assert fields[:3] == [f"rs{index}", "1", str(1000 + index)]
        assert math.isclose(float(fields[3]), min(frequency, 1.0 - frequency), abs_tol=1e-6)
        assert math.isclose(float(fields[4]), sum(pair_values) / len(pair_values), rel_tol=2e-5, abs_tol=1e-6)
        assert fields[5] == str(len(MARKERS) - 1)
        assert math.isclose(float(fields[6]), max(pair_values), abs_tol=1e-6)
        assert math.isclose(float(fields[7]), 1.0 + sum(pair_values), rel_tol=2e-5, abs_tol=1e-6)
    return output.with_suffix(".score.ld")


def check_pruning(gcta, work_dir, threads):
    prefix = Path(work_dir) / "input" / "ld"
    write_plink(prefix, MARKERS)
    output = run_gcta(
        gcta,
        work_dir,
        "prune",
        ["--bfile", prefix, "--ld-pruning", 0.5, "--ld-wind", 1000],
        threads=threads,
    )
    expected = expected_rsq()
    removed = {
        later
        for later in range(len(MARKERS))
        for earlier in range(later)
        if expected[later][earlier] > 0.5
    }
    expected_ids = [f"rs{index}" for index in range(len(MARKERS)) if index not in removed]
    actual_ids = output.with_suffix(".prune.in").read_text(encoding="utf-8").split()
    assert actual_ids == expected_ids
    return output.with_suffix(".prune.in")


def main():
    parser = argparse.ArgumentParser(description="GCTA legacy LD regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("score", "pruning"))
    args = parser.parse_args()
    if args.scenario == "score":
        single = check_score(args.gcta, Path(args.work) / "single", 1)
        multi = check_score(args.gcta, Path(args.work) / "multi", 4)
    else:
        single = check_pruning(args.gcta, Path(args.work) / "single", 1)
        multi = check_pruning(args.gcta, Path(args.work) / "multi", 4)
    assert_tables_close(single, multi)
    print(f"Passed legacy LD scenario: {args.scenario}")


if __name__ == "__main__":
    main()
