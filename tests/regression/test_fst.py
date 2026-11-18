#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink

MARKERS = [
    [0, 0, 0, 0, 2, 2, 2, 2],
]
GROUPS = ["A", "A", "A", "A", "B", "B", "B", "B"]


def expected_fst(first, second):
    sample_counts = [4.0, 4.0]
    n_bar = 4.0
    n_bar_r = 8.0
    n_c = 6.0
    p_bar = (first + second) / 2.0
    s_sq = (first - p_bar) ** 2 + (second - p_bar) ** 2
    h_bar = first * (1.0 - first) + second * (1.0 - second)
    d_buf = p_bar * (1.0 - p_bar) - 0.5 * s_sq
    a = (n_bar / n_c) * (s_sq - (d_buf - 0.25 * h_bar) / (n_bar - 1.0))
    b = (n_bar / (n_bar - 1.0)) * (d_buf - (2 * n_bar - 1.0) * h_bar / (4 * n_bar))
    c = 0.5 * h_bar
    return a / (a + b + c)


def main():
    parser = argparse.ArgumentParser(description="GCTA Fst regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", default="output", choices=("output",))
    args = parser.parse_args()
    outputs = []
    for threads in (1, 4):
        work = Path(args.work) / f"threads_{threads}"
        prefix = work / "input" / "populations"
        write_plink(prefix, MARKERS)
        subpopu = work / "subpopu.txt"
        subpopu.write_text(
            "".join(f"F1\tI{index}\t{group}\n" for index, group in enumerate(GROUPS)),
            encoding="utf-8",
        )
        output = run_gcta(
            args.gcta,
            work,
            "fst",
            ["--bfile", prefix, "--fst", "--sub-popu", subpopu],
            threads=threads,
        )
        rows = read_table(output.with_suffix(".fst"))
        assert len(rows) == 2
        assert rows[0][-1] == "Fst"
        for index, fields in enumerate(rows[1:]):
            first = MARKERS[index][0] / 2.0
            second = MARKERS[index][4] / 2.0
            assert math.isclose(float(fields[4]), first, abs_tol=1e-6)
            assert math.isclose(float(fields[5]), second, abs_tol=1e-6)
            assert math.isclose(float(fields[6]), expected_fst(first, second), rel_tol=2e-5, abs_tol=1e-6)
        outputs.append(output.with_suffix(".fst"))
    assert_tables_close(*outputs)
    print("Passed Fst frequency and variance oracles with 1/4-thread agreement")


if __name__ == "__main__":
    main()
