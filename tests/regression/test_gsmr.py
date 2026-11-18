#!/usr/bin/env python3

import argparse
import gzip
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink
from test_errors import expect_error
from test_mtcojo import summary_statistics as overlap_summary_statistics, write_ld_scores

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


def expected_pair(exposure, outcome, ratio):
    variances = [
        (outcome[index]["se"] / exposure[index]["b"]) ** 2
        + ratio**2 * (exposure[index]["se"] / exposure[index]["b"]) ** 2
        + 1e-6
        for index in range(len(exposure))
    ]
    precision = sum(1.0 / variance for variance in variances)
    se = math.sqrt(1.0 / precision)
    p = math.erfc(abs(ratio / se) / math.sqrt(2.0))
    return ratio, se, p


def expected_result(reverse=False, indices=None):
    ratio = 2.0 / 3.0 if reverse else 1.5
    statistics = summary_statistics(reverse)
    if indices is not None:
        statistics = [statistics[index] for index in indices]
    se = math.sqrt(1.0 / sum(1.0 / stat["variance"] for stat in statistics))
    p = math.erfc(abs(ratio / se) / math.sqrt(2.0))
    return ratio, se, p


def compress_summary(source, destination, flip_index=None):
    lines = Path(source).read_text(encoding="utf-8").splitlines()
    if flip_index is not None:
        fields = lines[flip_index + 1].split()
        fields[1], fields[2] = fields[2], fields[1]
        fields[3] = str(1.0 - float(fields[3]))
        fields[4] = str(-float(fields[4]))
        lines[flip_index + 1] = " ".join(fields)
    with gzip.open(destination, "wt", encoding="utf-8") as stream:
        stream.write("\n".join(lines) + "\n")


def run_output(gcta, work_dir, threads, compressed_flip=False):
    prefix, exposures, outcomes = make_fixture(work_dir)
    if compressed_flip:
        exposure = Path(exposures.read_text(encoding="utf-8").split()[1])
        outcome = Path(outcomes.read_text(encoding="utf-8").split()[1])
        exposure_gz = exposure.with_suffix(".ma.gz")
        outcome_gz = outcome.with_suffix(".ma.gz")
        compress_summary(exposure, exposure_gz)
        compress_summary(outcome, outcome_gz, flip_index=3)
        exposures.write_text(f"EXP {exposure_gz}\n", encoding="utf-8")
        outcomes.write_text(f"OUT {outcome_gz}\n", encoding="utf-8")
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
    assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
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


def write_custom_summary(path, statistics):
    Path(path).write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G 0.5 {stat['b']} {stat['se']} {stat['p']} 1000\n"
            for index, stat in enumerate(statistics)
        ),
        encoding="utf-8",
    )


def scaled_statistics(beta_scale, se_scale):
    base = summary_statistics(reverse=True)
    statistics = []
    for stat in base:
        b = stat["b"] * beta_scale
        se = stat["se"] * se_scale
        statistics.append({"b": b, "se": se, "p": math.erfc(abs(b / se) / math.sqrt(2.0))})
    return statistics


def run_multi_trait(gcta, work_dir, threads):
    work_dir = Path(work_dir)
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, hadamard_markers())
    exposure1 = scaled_statistics(1.0, 1.0)
    exposure2 = scaled_statistics(2.0, 1.0)
    outcome1 = scaled_statistics(1.5, 1.5)
    outcome2 = scaled_statistics(-1.0, 1.25)
    paths = {}
    for name, statistics in (("E1", exposure1), ("E2", exposure2), ("O1", outcome1), ("O2", outcome2)):
        path = work_dir / f"{name}.ma"
        write_custom_summary(path, statistics)
        paths[name] = path
    exposures = work_dir / "exposures.txt"
    outcomes = work_dir / "outcomes.txt"
    exposures.write_text(f"E1 {paths['E1']}\nE2 {paths['E2']}\n", encoding="utf-8")
    outcomes.write_text(f"O1 {paths['O1']}\nO2 {paths['O2']}\n", encoding="utf-8")
    output = run_gcta(
        gcta,
        work_dir,
        "multi",
        ["--bfile", prefix, "--gsmr-file", exposures, outcomes, "--gsmr-snp-min", 5, "--heidi-thresh", 0],
        threads=threads,
    )
    result = output.with_suffix(".gsmr")
    assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
    rows = read_table(result)
    assert len(rows) == 5
    expected = {
        ("E1", "O1"): (exposure1, outcome1, 1.5),
        ("E1", "O2"): (exposure1, outcome2, -1.0),
        ("E2", "O1"): (exposure2, outcome1, 0.75),
        ("E2", "O2"): (exposure2, outcome2, -0.5),
    }
    for row in rows[1:]:
        exposure, outcome, ratio = expected[tuple(row[:2])]
        expected_b, expected_se, expected_p = expected_pair(exposure, outcome, ratio)
        assert math.isclose(float(row[2]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[3]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[4]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
        assert row[5] == str(MARKER_COUNT)
    return result


def run_correlated(gcta, work_dir, threads):
    work_dir = Path(work_dir)
    markers = hadamard_markers()
    markers[-1] = list(markers[-2])
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, markers)
    exposure = work_dir / "exposure.ma"
    outcome = work_dir / "outcome.ma"
    write_summary(exposure, reverse=True)
    write_summary(outcome)
    exposures = work_dir / "exposures.txt"
    outcomes = work_dir / "outcomes.txt"
    exposures.write_text(f"EXP {exposure}\n", encoding="utf-8")
    outcomes.write_text(f"OUT {outcome}\n", encoding="utf-8")
    output = run_gcta(
        gcta,
        work_dir,
        "correlated",
        ["--bfile", prefix, "--gsmr-file", exposures, outcomes, "--gsmr-snp-min", 5, "--heidi-thresh", 0],
        threads=threads,
    )
    result = output.with_suffix(".gsmr")
    assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
    rows = read_table(result)
    expected_b, expected_se, expected_p = expected_result(indices=range(MARKER_COUNT - 1))
    assert math.isclose(float(rows[1][2]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
    assert math.isclose(float(rows[1][3]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
    assert math.isclose(float(rows[1][4]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
    assert rows[1][5] == str(MARKER_COUNT - 1)
    return result


def run_heidi(gcta, work_dir, threads):
    work_dir = Path(work_dir)
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, hadamard_markers())
    exposure_stats = summary_statistics(reverse=True)
    outcome_stats = summary_statistics()
    outlier = 5
    outcome_stats[outlier]["b"] = 10.0 * exposure_stats[outlier]["b"]
    outcome_stats[outlier]["se"] = 0.0001
    outcome_stats[outlier]["p"] = math.erfc(
        abs(outcome_stats[outlier]["b"] / outcome_stats[outlier]["se"]) / math.sqrt(2.0)
    )
    exposure = work_dir / "exposure.ma"
    outcome = work_dir / "outcome.ma"
    write_custom_summary(exposure, exposure_stats)
    write_custom_summary(outcome, outcome_stats)
    exposures = work_dir / "exposures.txt"
    outcomes = work_dir / "outcomes.txt"
    exposures.write_text(f"EXP {exposure}\n", encoding="utf-8")
    outcomes.write_text(f"OUT {outcome}\n", encoding="utf-8")
    output = run_gcta(
        gcta,
        work_dir,
        "heidi",
        ["--bfile", prefix, "--gsmr-file", exposures, outcomes, "--gsmr-snp-min", 10],
        threads=threads,
    )
    result = output.with_suffix(".gsmr")
    assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
    rows = read_table(result)
    retained = [index for index in range(MARKER_COUNT) if index != outlier]
    expected_b, expected_se, expected_p = expected_result(indices=retained)
    assert math.isclose(float(rows[1][2]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
    assert math.isclose(float(rows[1][3]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
    assert math.isclose(float(rows[1][4]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
    assert rows[1][5] == str(MARKER_COUNT - 1)
    pleio = Path(f"{output}.pleio_snps").read_text(encoding="utf-8")
    assert pleio == f"EXP OUT rs{outlier}\n"
    return result


def write_overlap_summary(path, trait):
    Path(path).write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G 0.5 {stat['b']} {stat['se']} {stat['p']} 1000\n"
            for index, stat in enumerate(overlap_summary_statistics(trait))
        ),
        encoding="utf-8",
    )


def make_overlap_fixture(work_dir):
    work_dir = Path(work_dir)
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, hadamard_markers())
    exposure = work_dir / "exposure.ma"
    outcome = work_dir / "outcome.ma"
    write_overlap_summary(exposure, "X")
    write_overlap_summary(outcome, "Y")
    exposures = work_dir / "exposures.txt"
    outcomes = work_dir / "outcomes.txt"
    exposures.write_text(f"EXP {exposure}\n", encoding="utf-8")
    outcomes.write_text(f"OUT {outcome}\n", encoding="utf-8")
    reference, weights = write_ld_scores(work_dir)
    return prefix, exposures, outcomes, reference, weights


def expected_overlap_result():
    exposure = overlap_summary_statistics("X")[:2]
    outcome = overlap_summary_statistics("Y")[:2]
    ratios = [outcome[index]["b"] / exposure[index]["b"] for index in range(2)]
    ordered = sorted(ratios)
    median = (ordered[0] + ordered[1]) / 2.0
    variances = [
        (outcome[index]["se"] / exposure[index]["b"]) ** 2
        + median**2 * (exposure[index]["se"] / exposure[index]["b"]) ** 2
        + 1e-6
        for index in range(2)
    ]
    precision = sum(1.0 / variance for variance in variances)
    estimate = sum(ratio / variance for ratio, variance in zip(ratios, variances)) / precision
    se = math.sqrt(1.0 / precision)
    p = math.erfc(abs(estimate / se) / math.sqrt(2.0))
    return estimate, se, p


def run_overlap(gcta, work_dir, algorithm, threads):
    prefix, exposures, outcomes, reference, weights = make_overlap_fixture(work_dir)
    output = run_gcta(
        gcta,
        work_dir,
        "overlap",
        [
            "--bfile",
            prefix,
            "--gsmr-file",
            exposures,
            outcomes,
            "--gsmr-so",
            algorithm,
            "--ref-ld-chr",
            reference,
            "--w-ld-chr",
            weights,
            "--gsmr-snp-min",
            2,
            "--heidi-thresh",
            0,
        ],
        threads=threads,
    )
    result = output.with_suffix(".gsmr")
    assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
    rows = read_table(result)
    assert len(rows) == 2
    expected_b, expected_se, expected_p = expected_overlap_result()
    assert math.isclose(float(rows[1][2]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
    assert math.isclose(float(rows[1][3]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
    assert math.isclose(float(rows[1][4]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
    assert rows[1][5] == "2"
    log = output.with_suffix(".log").read_text(encoding="utf-8")
    if algorithm == 0:
        assert "LD score regression analysis completed." in log
    else:
        exposure = overlap_summary_statistics("X")
        outcome = overlap_summary_statistics("Y")
        null_indices = [
            index for index in range(MARKER_COUNT)
            if exposure[index]["p"] >= 0.01 and outcome[index]["p"] >= 0.01
        ]
        assert len(null_indices) == MARKER_COUNT - 2
        x_values = [exposure[index]["z"] for index in null_indices]
        y_values = [outcome[index]["z"] for index in null_indices]
        x_mean = sum(x_values) / len(x_values)
        y_mean = sum(y_values) / len(y_values)
        x_centered = [value - x_mean for value in x_values]
        y_centered = [value - y_mean for value in y_values]
        expected_correlation = sum(x * y for x, y in zip(x_centered, y_centered)) / math.sqrt(
            sum(value * value for value in x_centered)
            * sum(value * value for value in y_centered)
        )
        lines = log.splitlines()
        correlation_line = lines[lines.index("Correlation:") + 1]
        observed_correlation = float(correlation_line.split()[1])
        assert math.isclose(observed_correlation, expected_correlation, rel_tol=2e-5, abs_tol=1e-6)
    return result


def expect_invalid_overlap(gcta, work_dir):
    work_dir = Path(work_dir)
    prefix, exposures, outcomes = make_fixture(work_dir)
    output = work_dir / "invalid"
    expect_error(
        gcta,
        work_dir,
        "invalid",
        [
            "--bfile",
            prefix,
            "--gsmr-file",
            exposures,
            outcomes,
            "--gsmr-so",
            2,
            "--out",
            output,
        ],
    )
    console = (work_dir / "invalid.console").read_text(encoding="utf-8")
    assert "--gsmr-so should be 0 (LD score regression) or 1" in console


def main():
    parser = argparse.ArgumentParser(description="GCTA GSMR regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument(
        "--scenario",
        default="output",
        choices=(
            "output",
            "compressed_flip",
            "multi_trait",
            "correlated",
            "heidi",
            "overlap_ldsc",
            "overlap_correlation",
            "invalid_overlap",
        ),
    )
    args = parser.parse_args()
    if args.scenario == "invalid_overlap":
        expect_invalid_overlap(args.gcta, args.work)
    elif args.scenario in ("multi_trait", "correlated", "heidi"):
        runner = {"multi_trait": run_multi_trait, "correlated": run_correlated, "heidi": run_heidi}[args.scenario]
        single = runner(args.gcta, Path(args.work) / "single", 1)
        multi = runner(args.gcta, Path(args.work) / "multi", 4)
        assert_tables_close(single, multi)
    elif args.scenario.startswith("overlap_"):
        algorithm = 0 if args.scenario == "overlap_ldsc" else 1
        single = run_overlap(args.gcta, Path(args.work) / "single", algorithm, 1)
        multi = run_overlap(args.gcta, Path(args.work) / "multi", algorithm, 4)
        assert_tables_close(single, multi)
    else:
        single = run_output(args.gcta, Path(args.work) / "single", 1, args.scenario == "compressed_flip")
        multi = run_output(args.gcta, Path(args.work) / "multi", 4, args.scenario == "compressed_flip")
        assert_tables_close(single, multi)
    print(f"Passed GSMR regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
