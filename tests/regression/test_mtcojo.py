#!/usr/bin/env python3

import argparse
import gzip
import math
import os
from pathlib import Path
import subprocess

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink

SAMPLE_COUNT = 16
MARKER_COUNT = 12
LD_SCORES = (1, 2, 4, 7, 11, 16)
BXY = 0.4


def orthogonal_markers(sample_count=SAMPLE_COUNT, marker_count=MARKER_COUNT):
    matrix = [[1]]
    while len(matrix) < sample_count:
        size = len(matrix)
        matrix = (
            [matrix[row] + matrix[row] for row in range(size)]
            + [matrix[row] + [-value for value in matrix[row]] for row in range(size)]
        )
    return [
        [0 if matrix[sample][marker] < 0 else 2 for sample in range(sample_count)]
        for marker in range(1, marker_count + 1)
    ]


def summary_statistics(trait):
    statistics = []
    for index in range(MARKER_COUNT):
        pair = index // 2
        ld_score = LD_SCORES[pair]
        z_magnitude = 6.0 if pair == 0 else math.sqrt(1.0 + (0.08 if trait == "Y" else 0.1) * ld_score)
        z_value = z_magnitude if trait == "Y" or index % 2 == 0 else -z_magnitude
        se = 0.04 + 0.001 * index if trait == "Y" else 0.05 + 0.002 * index
        b = z_value * se
        p = math.erfc(abs(z_value) / math.sqrt(2.0))
        statistics.append({"b": b, "se": se, "p": p, "z": z_value})
    return statistics


def write_summary(path, trait):
    Path(path).write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G 0.5 {stat['b']} {stat['se']} {stat['p']} 1000\n"
            for index, stat in enumerate(summary_statistics(trait))
        ),
        encoding="utf-8",
    )


def write_ld_scores(work_dir, scores=LD_SCORES, marker_count=MARKER_COUNT):
    reference = Path(work_dir) / "reference_ld"
    weights = Path(work_dir) / "regression_weights"
    for directory in (reference, weights):
        directory.mkdir(parents=True, exist_ok=True)
        for chromosome in range(1, 23):
            (directory / f"{chromosome}.l2.M_5_50").write_text("100\n", encoding="utf-8")
            weighted = directory == weights
            lines = ["CHR SNP BP CM MAF L2" if weighted else "CHR SNP BP L2"]
            if chromosome == 1:
                for index in range(marker_count):
                    ld_score = scores[index // 2]
                    if weighted:
                        lines.append(f"1 rs{index} {1000 + index} 0 0.5 {1.0 + 0.1 * ld_score}")
                    else:
                        lines.append(f"1 rs{index} {1000 + index} {ld_score}")
            (directory / f"{chromosome}.l2.ldscore").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return reference, weights


def make_fixture(work_dir):
    work_dir = Path(work_dir)
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, orthogonal_markers())
    target = work_dir / "target.ma"
    covariate = work_dir / "covariate.ma"
    write_summary(target, "Y")
    write_summary(covariate, "X")
    trait_list = work_dir / "traits.txt"
    trait_list.write_text(f"Y {target}\nX {covariate}\n", encoding="utf-8")
    bxy = work_dir / "bxy.txt"
    bxy.write_text(f"X {BXY}\n", encoding="utf-8")
    reference, weights = write_ld_scores(work_dir)
    return prefix, trait_list, bxy, reference, weights


def mtcojo_arguments(prefix, trait_list, bxy, reference, weights):
    return [
        "--bfile",
        prefix,
        "--mtcojo-file",
        trait_list,
        "--mtcojo-bxy",
        bxy,
        "--ref-ld-chr",
        reference,
        "--w-ld-chr",
        weights,
    ]


def check_zero_cross_intercept(output):
    lines = output.with_suffix(".log").read_text(encoding="utf-8").splitlines()
    start = lines.index("Intercept:")
    values = [[float(value) for value in line.split()] for line in lines[start + 1:start + 3]]
    assert values[0][0] == 1.0 and values[1][1] == 1.0
    assert values[0][1] == 0.0 and values[1][0] == 0.0


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
    prefix, trait_list, bxy, reference, weights = make_fixture(work_dir)
    if compressed_flip:
        target = Path(trait_list.read_text(encoding="utf-8").splitlines()[0].split()[1])
        covariate = Path(trait_list.read_text(encoding="utf-8").splitlines()[1].split()[1])
        target_gz = target.with_suffix(".ma.gz")
        covariate_gz = covariate.with_suffix(".ma.gz")
        compress_summary(target, target_gz)
        compress_summary(covariate, covariate_gz, flip_index=4)
        trait_list.write_text(f"Y {target_gz}\nX {covariate_gz}\n", encoding="utf-8")
    output = run_gcta(
        gcta,
        work_dir,
        "mtcojo",
        mtcojo_arguments(prefix, trait_list, bxy, reference, weights),
        threads=threads,
    )
    check_zero_cross_intercept(output)
    result = Path(f"{output}.mtcojo.cma")
    rows = read_table(result)
    assert rows[0] == ["SNP", "A1", "A2", "freq", "b", "se", "p", "N", "bC", "bC_se", "bC_pval"]
    assert len(rows) == MARKER_COUNT + 1
    target = summary_statistics("Y")
    covariate = summary_statistics("X")
    for index, row in enumerate(rows[1:]):
        assert row[:2] == [f"rs{index}", "A"]
        assert row[2:4] == ["G", "0.5"]
        assert math.isclose(float(row[4]), target[index]["b"], rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[5]), target[index]["se"], rel_tol=2e-5, abs_tol=1e-6)
        expected_b = target[index]["b"] - BXY * covariate[index]["b"]
        expected_se = math.sqrt(target[index]["se"] ** 2 + BXY**2 * covariate[index]["se"] ** 2)
        expected_p = math.erfc(abs(expected_b / expected_se) / math.sqrt(2.0))
        assert math.isclose(float(row[8]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[9]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[10]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
    return result


INTERNAL_COUNT = 40
INTERNAL_LD_SCORES = tuple(1 + (index * 7 % 20) for index in range(INTERNAL_COUNT // 2))


def internal_statistics(trait):
    statistics = []
    for index in range(INTERNAL_COUNT):
        ld_score = INTERNAL_LD_SCORES[index // 2]
        z_magnitude = 5.5 if index < 10 else math.sqrt(1.0 + (0.08 if trait == "Y" else 0.1) * ld_score)
        z_value = z_magnitude if trait == "Y" or index % 2 == 0 else -z_magnitude
        se = 0.04 + 0.001 * (index % 12) if trait == "Y" else 0.05 + 0.002 * (index % 12)
        b = z_value * se
        statistics.append({"b": b, "se": se, "p": math.erfc(abs(z_value) / math.sqrt(2.0)), "z": z_value})
    return statistics


def write_internal_summary(path, trait):
    Path(path).write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G 0.5 {stat['b']} {stat['se']} {stat['p']} 1000\n"
            for index, stat in enumerate(internal_statistics(trait))
        ),
        encoding="utf-8",
    )


def expected_internal_bxy():
    exposure = internal_statistics("X")[:10]
    outcome = internal_statistics("Y")[:10]
    ratios = [outcome[index]["b"] / exposure[index]["b"] for index in range(10)]
    ordered = sorted(ratios)
    median = (ordered[4] + ordered[5]) / 2.0
    variances = [
        (outcome[index]["se"] / exposure[index]["b"]) ** 2
        + median**2 * (exposure[index]["se"] / exposure[index]["b"]) ** 2
        + 1e-6
        for index in range(10)
    ]
    precision = sum(1.0 / variance for variance in variances)
    return sum(ratio / variance for ratio, variance in zip(ratios, variances)) / precision


def run_internal_bxy(gcta, work_dir, threads):
    work_dir = Path(work_dir)
    prefix = work_dir / "input" / "reference"
    write_plink(prefix, orthogonal_markers(64, INTERNAL_COUNT))
    target = work_dir / "target.ma"
    covariate = work_dir / "covariate.ma"
    write_internal_summary(target, "Y")
    write_internal_summary(covariate, "X")
    trait_list = work_dir / "traits.txt"
    trait_list.write_text(f"Y {target}\nX {covariate}\n", encoding="utf-8")
    reference, weights = write_ld_scores(work_dir, INTERNAL_LD_SCORES, INTERNAL_COUNT)
    output = run_gcta(
        gcta,
        work_dir,
        "internal",
        ["--bfile", prefix, "--mtcojo-file", trait_list, "--ref-ld-chr", reference, "--w-ld-chr", weights, "--heidi-thresh", 0],
        threads=threads,
    )
    check_zero_cross_intercept(output)
    result = Path(f"{output}.mtcojo.cma")
    rows = read_table(result)
    assert len(rows) == INTERNAL_COUNT + 1
    expected_bxy = expected_internal_bxy()
    target_stats = internal_statistics("Y")
    covariate_stats = internal_statistics("X")
    for index, row in enumerate(rows[1:]):
        expected_b = target_stats[index]["b"] - expected_bxy * covariate_stats[index]["b"]
        expected_se = math.sqrt(target_stats[index]["se"] ** 2 + expected_bxy**2 * covariate_stats[index]["se"] ** 2)
        expected_p = math.erfc(abs(expected_b / expected_se) / math.sqrt(2.0))
        assert math.isclose(float(row[8]), expected_b, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[9]), expected_se, rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(row[10]), expected_p, rel_tol=2e-5, abs_tol=1e-100)
    return result


def expect_missing_ld(gcta, work_dir):
    work_dir = Path(work_dir)
    prefix, trait_list, bxy, _, _ = make_fixture(work_dir)
    output = work_dir / "missing_ld"
    arguments = mtcojo_arguments(
        prefix,
        trait_list,
        bxy,
        work_dir / "absent_reference",
        work_dir / "absent_weights",
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    command = [
        str(Path(gcta).resolve()),
        *(str(argument.resolve() if isinstance(argument, Path) else argument) for argument in arguments),
        "--out",
        str(output),
        "--thread-num",
        "1",
    ]
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
    (work_dir / "missing_ld.console").write_text(result.stdout, encoding="utf-8")
    assert result.returncode != 0
    assert "cannot open the file" in result.stdout
    assert "1.l2.M_5_50" in result.stdout


def main():
    parser = argparse.ArgumentParser(description="GCTA mtCOJO regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("output", "compressed_flip", "internal_bxy", "missing_ld"))
    args = parser.parse_args()
    if args.scenario == "missing_ld":
        expect_missing_ld(args.gcta, args.work)
    elif args.scenario == "internal_bxy":
        single = run_internal_bxy(args.gcta, Path(args.work) / "single", 1)
        multi = run_internal_bxy(args.gcta, Path(args.work) / "multi", 4)
        assert_tables_close(single, multi)
    else:
        single = run_output(args.gcta, Path(args.work) / "single", 1, args.scenario == "compressed_flip")
        multi = run_output(args.gcta, Path(args.work) / "multi", 4, args.scenario == "compressed_flip")
        assert_tables_close(single, multi)
    print(f"Passed mtCOJO regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
