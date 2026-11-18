#!/usr/bin/env python3

import argparse
import gzip
import math
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, assert_tables_close, read_floats, read_table, run_gcta
from test_errors import expect_error

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


def run_erm_gz(gcta, work_dir, threads):
    expression = Path(work_dir) / "input" / "expression.txt"
    write_expression(expression)
    output = run_gcta(gcta, work_dir, "erm", ["--efile", expression, "--make-erm-gz"], threads=threads)
    result = Path(f"{output}.E.grm.gz")
    expected, expected_denominators = expected_erm(0)
    with gzip.open(result, "rt", encoding="utf-8") as stream:
        rows = [line.split() for line in stream.read().splitlines()]
    assert len(rows) == len(expected)
    index = 0
    for row in range(len(PROBE_DATA)):
        for column in range(row + 1):
            assert rows[index][:2] == [str(row + 1), str(column + 1)]
            assert math.isclose(float(rows[index][2]), expected_denominators[index], rel_tol=2e-5, abs_tol=1e-6)
            assert math.isclose(float(rows[index][3]), expected[index], rel_tol=2e-5, abs_tol=1e-6)
            index += 1
    ids = Path(f"{output}.E.grm.id").read_text(encoding="utf-8").splitlines()
    assert ids == [f"F{index}\tI{index}" for index in range(len(PROBE_DATA))]
    return result


def expect_duplicate_probe(gcta, work_dir):
    expression = Path(work_dir) / "input" / "duplicate.txt"
    expression.parent.mkdir(parents=True, exist_ok=True)
    expression.write_text("FID IID P0 P0\nF0 I0 1 2\n", encoding="utf-8")
    output = Path(work_dir) / "duplicate"
    expect_error(gcta, Path(work_dir), "duplicate", ["--efile", expression, "--make-erm", "--out", output])
    console = (Path(work_dir) / "duplicate.console").read_text(encoding="utf-8")
    assert "duplicated probe names found" in console


ECOJO_DATA = [
    [1.0, 0.5],
    [0.0, 1.0],
    [-1.0, -0.5],
    [0.0, -1.0],
]
ECOJO_Z = [30.0, 30.0]
ECOJO_N = [100.0, 100.0]
ECOJO_LAMBDA = 0.5


def ecojo_correlation():
    centered = []
    for probe in range(2):
        values = [row[probe] for row in ECOJO_DATA]
        mean = sum(values) / len(values)
        centered.append([value - mean for value in values])
    numerator = sum(left * right for left, right in zip(centered[0], centered[1]))
    denominator = math.sqrt(
        sum(value * value for value in centered[0])
        * sum(value * value for value in centered[1])
    )
    return numerator / denominator


def write_ecojo_fixture(work_dir):
    work_dir = Path(work_dir)
    expression = work_dir / "input" / "ecojo.txt"
    expression.parent.mkdir(parents=True, exist_ok=True)
    expression.write_text(
        "FID IID P0 P1\n"
        + "".join(
            f"F{index} I{index} {row[0]} {row[1]}\n" for index, row in enumerate(ECOJO_DATA)
        ),
        encoding="utf-8",
    )
    summary = work_dir / "ecojo.ma"
    summary.write_text(
        "Probe z N\n"
        + "".join(f"P{index} {ECOJO_Z[index]} {ECOJO_N[index]}\n" for index in range(2)),
        encoding="utf-8",
    )
    correlation = work_dir / "ecojo.cor"
    r_value = ecojo_correlation()
    correlation.write_text(f"P0 P1\n1 {r_value}\n{r_value} 1\n", encoding="utf-8")
    return expression, summary, correlation


def solve_2x2(matrix, vector):
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    return [
        (matrix[1][1] * vector[0] - matrix[0][1] * vector[1]) / determinant,
        (matrix[0][0] * vector[1] - matrix[1][0] * vector[0]) / determinant,
    ]


def ecojo_marginals():
    scales = [math.sqrt(ECOJO_Z[index] ** 2 + ECOJO_N[index]) for index in range(2)]
    b = [ECOJO_Z[index] / scales[index] for index in range(2)]
    se = [1.0 / scales[index] for index in range(2)]
    return b, se


def check_ecojo_blup(result, expected_b, expected_se, expected_blup):
    rows = read_table(result)
    assert rows[0] == ["Probe", "b", "se", "z", "n", "b_blup"]
    assert len(rows) == 3
    for index, row in enumerate(rows[1:]):
        assert row[0] == f"P{index}"
        assert math.isclose(float(row[1]), expected_b[index], rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[2]), expected_se[index], rel_tol=2e-5, abs_tol=1e-6)
        assert row[3] == str(int(ECOJO_Z[index]))
        assert row[4] == str(int(ECOJO_N[index]))
        assert math.isclose(float(row[5]), expected_blup[index], rel_tol=2e-5, abs_tol=1e-6)


def run_ecojo_select(gcta, work_dir, threads):
    expression, summary, _ = write_ecojo_fixture(work_dir)
    output = run_gcta(
        gcta,
        work_dir,
        "select",
        ["--efile", expression, "--ecojo", summary, "--ecojo-slct"],
        threads=threads,
    )
    result = Path(f"{output}.slct.ecojo")
    rows = read_table(result)
    assert len(rows) == 3
    b, se = ecojo_marginals()
    r_value = ecojo_correlation()
    correlation = [[1.0, r_value], [r_value, 1.0]]
    expected_bj = solve_2x2(correlation, b)
    inverse_diagonal = 1.0 / (1.0 - r_value * r_value)
    for index, row in enumerate(rows[1:]):
        assert row[0] == f"P{index}"
        assert math.isclose(float(row[1]), b[index], rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[2]), se[index], rel_tol=2e-5, abs_tol=1e-6)
        expected_sej = math.sqrt(inverse_diagonal / ECOJO_N[index])
        expected_p = math.erfc(abs(expected_bj[index] / expected_sej) / math.sqrt(2.0))
        assert math.isclose(float(row[5]), expected_bj[index], rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[6]), expected_sej, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[7]), expected_bj[index] / expected_sej, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[8]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
    return result


def run_ecojo_blup(gcta, work_dir, threads):
    expression, summary, correlation = write_ecojo_fixture(work_dir)
    outputs = []
    for name, option, data in (
        ("blup_file", "--efile", expression),
        ("blup_cor", "--e-cor", correlation),
    ):
        output = run_gcta(
            gcta,
            work_dir,
            name,
            [option, data, "--ecojo", summary, "--ecojo-blup", ECOJO_LAMBDA],
            threads=threads,
        )
        outputs.append(Path(f"{output}.blup.ecojo"))
    b, se = ecojo_marginals()
    r_value = ecojo_correlation()
    diagonal = 1.0 + 2.0 * (1.0 / ECOJO_LAMBDA - 1.0) / (sum(ECOJO_N) / 2.0)
    expected_blup = solve_2x2([[diagonal, r_value], [r_value, diagonal]], b)
    for result in outputs:
        check_ecojo_blup(result, b, se, expected_blup)
    assert_tables_close(*outputs)
    return outputs[0]


def run_ecojo_no_signal(gcta, work_dir, threads):
    expression, summary, _ = write_ecojo_fixture(work_dir)
    summary.write_text("Probe z N\nP0 0 100\nP1 0 100\n", encoding="utf-8")
    output = run_gcta(
        gcta,
        work_dir,
        "no_signal",
        ["--efile", expression, "--ecojo", summary, "--ecojo-slct"],
        threads=threads,
    )
    assert not Path(f"{output}.slct.ecojo").exists()
    console = Path(f"{output}.console").read_text(encoding="utf-8")
    assert "No probe has been selected." in console


def main():
    parser = argparse.ArgumentParser(description="GCTA expression workflow regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument(
        "--scenario",
        required=True,
        choices=(
            "erm_std",
            "erm_center",
            "erm_individual",
            "erm_gz",
            "duplicate_probe",
            "ecojo_select",
            "ecojo_blup",
            "ecojo_no_signal",
        ),
    )
    args = parser.parse_args()
    if args.scenario == "duplicate_probe":
        expect_duplicate_probe(args.gcta, args.work)
    elif args.scenario == "erm_gz":
        single = run_erm_gz(args.gcta, Path(args.work) / "single", 1)
        multi = run_erm_gz(args.gcta, Path(args.work) / "multi", 4)
        with gzip.open(single, "rt", encoding="utf-8") as stream:
            single_text = stream.read()
        with gzip.open(multi, "rt", encoding="utf-8") as stream:
            multi_text = stream.read()
        assert single_text == multi_text
    elif args.scenario.startswith("erm_"):
        single = run_erm(args.gcta, Path(args.work) / "single", args.scenario, 1)
        multi = run_erm(args.gcta, Path(args.work) / "multi", args.scenario, 4)
        assert_grm_equal(single, multi)
    elif args.scenario == "ecojo_no_signal":
        run_ecojo_no_signal(args.gcta, Path(args.work) / "single", 1)
        run_ecojo_no_signal(args.gcta, Path(args.work) / "multi", 4)
    else:
        runner = run_ecojo_select if args.scenario == "ecojo_select" else run_ecojo_blup
        single = runner(args.gcta, Path(args.work) / "single", 1)
        multi = runner(args.gcta, Path(args.work) / "multi", 4)
        assert_tables_close(single, multi)
    print(f"Passed expression regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
