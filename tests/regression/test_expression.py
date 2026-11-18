#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, run_gcta

PROBE_DATA = [
    [1.0, 2.0, None, 4.0],
    [2.0, 0.0, 3.0, None],
    [3.0, 1.0, 2.0, 2.0],
    [4.0, 3.0, 1.0, 0.0],
    [5.0, 2.0, 0.0, 1.0],
]


def write_expression(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "FID IID P0 P1 P2 P3\n"
        + "".join(
            f"F{index} I{index} " + " ".join("-9" if value is None else str(value) for value in row) + "\n"
            for index, row in enumerate(PROBE_DATA)
        ),
        encoding="utf-8",
    )


def expected_erm(algorithm):
    sample_count = len(PROBE_DATA)
    probe_count = len(PROBE_DATA[0])
    transformed = [[0.0] * probe_count for _ in range(sample_count)]
    variances = [0.0] * probe_count
    if algorithm < 2:
        for probe in range(probe_count):
            observed = [row[probe] for row in PROBE_DATA if row[probe] is not None]
            mean = sum(observed) / len(observed)
            centered = [(row[probe] - mean) if row[probe] is not None else 0.0 for row in PROBE_DATA]
            scale = math.sqrt(sum(value * value for value in centered) / (len(observed) - 1))
            for sample in range(sample_count):
                transformed[sample][probe] = centered[sample] / scale if algorithm == 0 else centered[sample]
            variances[probe] = sum(value * value for value in centered) / (len(observed) - 1)
    else:
        for sample, row in enumerate(PROBE_DATA):
            observed = [value for value in row if value is not None]
            mean = sum(observed) / len(observed)
            transformed[sample] = [(value - mean) if value is not None else 0.0 for value in row]

    raw = []
    denominators = []
    diagonal = []
    for row in range(sample_count):
        for column in range(row + 1):
            common = [
                probe for probe in range(probe_count)
                if PROBE_DATA[row][probe] is not None and PROBE_DATA[column][probe] is not None
            ]
            numerator = sum(transformed[row][probe] * transformed[column][probe] for probe in common)
            if algorithm == 0:
                denominator = float(len(common))
            elif algorithm == 1:
                denominator = sum(variances[probe] for probe in common)
            else:
                row_ssx = sum(transformed[row][probe] ** 2 for probe in common)
                column_ssx = sum(transformed[column][probe] ** 2 for probe in common)
                denominator = math.sqrt(row_ssx * column_ssx)
            raw.append(numerator / denominator if denominator > 0.0 else 0.0)
            denominators.append(denominator)
            if row == column:
                diagonal.append(raw[-1])
    scale = sum(diagonal) / sample_count
    return [value / scale for value in raw], denominators


def run_erm(gcta, work_dir, scenario, threads):
    expression = Path(work_dir) / "input" / "expression.txt"
    write_expression(expression)
    algorithm = {"erm_std": 0, "erm_center": 1, "erm_individual": 2}[scenario]
    option = ["--make-erm"] if algorithm == 0 else ["--make-erm-alg", str(algorithm + 1)]
    output = run_gcta(gcta, work_dir, "erm", ["--efile", expression, *option], threads=threads)
    prefix = Path(f"{output}.E")
    expected, expected_denominators = expected_erm(algorithm)
    assert_floats_close(read_floats(Path(f"{prefix}.grm.bin")), expected, prefix)
    assert_floats_close(read_floats(Path(f"{prefix}.grm.N.bin")), expected_denominators, prefix)
    ids = Path(f"{prefix}.grm.id").read_text(encoding="utf-8").splitlines()
    assert ids == [f"F{index}\tI{index}" for index in range(len(PROBE_DATA))]
    return prefix


def main():
    parser = argparse.ArgumentParser(description="GCTA expression relationship matrix regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("erm_std", "erm_center", "erm_individual"))
    args = parser.parse_args()
    single = run_erm(args.gcta, Path(args.work) / "single", args.scenario, 1)
    multi = run_erm(args.gcta, Path(args.work) / "multi", args.scenario, 4)
    assert_grm_equal(single, multi)
    print(f"Passed expression regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
