#!/usr/bin/env python3

import argparse
import gzip
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink

SAMPLE_COUNT = 16
MARKER_COUNT = 12


def hadamard_markers():
    matrix = [[1]]
    while len(matrix) < SAMPLE_COUNT:
        size = len(matrix)
        matrix = (
            [matrix[row] + matrix[row] for row in range(size)]
            + [matrix[row] + [-value for value in matrix[row]] for row in range(size)]
        )
    return [
        [0 if matrix[sample][marker] < 0 else 2 for sample in range(SAMPLE_COUNT)]
        for marker in range(1, MARKER_COUNT + 1)
    ]


def summary_statistics(reverse=False):
    statistics = []
    for index in range(MARKER_COUNT):
        exposure_b = 0.2 + 0.01 * index
        exposure_se = 0.02 + 0.001 * index
        outcome_b = 1.5 * exposure_b
        outcome_se = 0.03 + 0.0015 * index
        if reverse:
            b, se = exposure_b, exposure_se
            x, x_se, y_se = outcome_b, outcome_se, exposure_se
            ratio = 2.0 / 3.0
        else:
            b, se = outcome_b, outcome_se
            x, x_se, y_se = exposure_b, exposure_se, outcome_se
            ratio = 1.5
        p = math.erfc(abs(b / se) / math.sqrt(2.0))
        variance = (y_se / x) ** 2 + ratio**2 * (x_se / x) ** 2 + 1e-6
        statistics.append({"b": b, "se": se, "p": p, "variance": variance})
    return statistics


def write_summary(path, reverse=False):
    statistics = summary_statistics(reverse)
    Path(path).write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G 0.5 {stat['b']} {stat['se']} {stat['p']} 1000\n"
            for index, stat in enumerate(statistics)
        ),
        encoding="utf-8",
    )


def make_fixture(work_dir):
    work_dir = Path(work_dir)
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, hadamard_markers())
    exposure = work_dir / "exposure.ma"
    outcome = work_dir / "outcome.ma"
    write_summary(exposure, reverse=True)
    write_summary(outcome)
    exposure_list = work_dir / "exposures.txt"
    outcome_list = work_dir / "outcomes.txt"
    exposure_list.write_text(f"EXP {exposure}\n", encoding="utf-8")
    outcome_list.write_text(f"OUT {outcome}\n", encoding="utf-8")
    return prefix, exposure_list, outcome_list


def expected_result(reverse=False):
    ratio = 2.0 / 3.0 if reverse else 1.5
    statistics = summary_statistics(reverse)
    se = math.sqrt(1.0 / sum(1.0 / stat["variance"] for stat in statistics))
    p = math.erfc(abs(ratio / se) / math.sqrt(2.0))
    return ratio, se, p


def run_output(gcta, work_dir, threads):
    prefix, exposures, outcomes = make_fixture(work_dir)
    output = run_gcta(
        gcta,
        work_dir,
        "gsmr",
        [
            "--bfile",
            prefix,
            "--gsmr-file",
            exposures,
            outcomes,
            "--gsmr-direction",
            2,
            "--gsmr-snp-min",
            5,
            "--heidi-thresh",
            0,
            "--effect-plot",
        ],
        threads=threads,
    )
    result = output.with_suffix(".gsmr")
    rows = read_table(result)
    assert len(rows) == 3
    assert [row[:2] for row in rows[1:]] == [["EXP", "OUT"], ["OUT", "EXP"]]
    for row, reverse in zip(rows[1:], (False, True)):
        expected_b, expected_se, expected_p = expected_result(reverse)
        assert math.isclose(float(row[2]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[3]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[4]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
        assert row[5] == str(MARKER_COUNT)

    with gzip.open(f"{output}.eff_plot.gz", "rt", encoding="utf-8") as stream:
        effect_plot = stream.read().splitlines()
    assert effect_plot[:4] == ["#trait_begin", "EXP", "OUT", "#trait_end"]
    effect_begin = effect_plot.index("#effect_begin")
    effect_end = effect_plot.index("#effect_end")
    assert effect_end - effect_begin - 1 == MARKER_COUNT
    exposure_stats = summary_statistics(reverse=True)
    outcome_stats = summary_statistics()
    for line in effect_plot[effect_begin + 1:effect_end]:
        fields = line.split()
        index = int(fields[0][2:])
        assert fields[1:4] == ["A", "G", "0.5"]
        observed = [float(value) for value in fields[4:]]
        expected = [
            exposure_stats[index]["b"],
            exposure_stats[index]["se"],
            outcome_stats[index]["b"],
            outcome_stats[index]["se"],
        ]
        assert all(math.isclose(value, reference, abs_tol=1e-12) for value, reference in zip(observed, expected))
    return result



def main():
    parser = argparse.ArgumentParser(description="GCTA GSMR regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", default="output", choices=("output",))
    args = parser.parse_args()
    single = run_output(args.gcta, Path(args.work) / "single", 1)
    multi = run_output(args.gcta, Path(args.work) / "multi", 4)
    assert_tables_close(single, multi)
    print(f"Passed GSMR regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
