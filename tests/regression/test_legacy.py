#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, run_gcta
from plink import synthetic_markers, write_plink
from test_errors import expect_error


def expected_legacy_grm(markers, algorithm):
    assert algorithm in (0, 1)
    sample_count = len(markers[0])
    element_count = sample_count * (sample_count + 1) // 2
    sums = [0.0] * element_count
    counts = [0] * element_count
    variance_sum = 0.0
    marker_count = 0

    for marker in markers:
        observed = [value for value in marker if value is not None]
        if not observed:
            continue
        frequency = sum(observed) / (2.0 * len(observed))
        variance = 2.0 * frequency * (1.0 - frequency)
        variance_sum += variance
        marker_count += 1
        if algorithm == 0:
            scale = math.sqrt(1.0 / variance) if variance >= 1e-50 else 0.0
            values = [
                (value - 2.0 * frequency) * scale if value is not None else 0.0
                for value in marker
            ]
        else:
            values = [value - 2.0 * frequency if value is not None else 0.0 for value in marker]
        for row in range(sample_count):
            offset = row * (row + 1) // 2
            for column in range(row + 1):
                index = offset + column
                sums[index] += values[row] * values[column]
                if marker[row] is not None and marker[column] is not None:
                    counts[index] += 1

    assert marker_count > 0
    weight = 1.0 if algorithm == 0 else marker_count / variance_sum
    grm = [sums[index] / count * weight if count else 0.0 for index, count in enumerate(counts)]
    if algorithm == 1:
        mean_variance = variance_sum / marker_count
        counts = [count * mean_variance for count in counts]
    return grm, [float(count) for count in counts]


def make_fixture(work_dir):
    markers = synthetic_markers()
    markers[1][0] = None
    markers[2][1] = None
    markers[3][2] = None
    prefix = Path(work_dir) / "input" / "synthetic"
    sample_ids = write_plink(prefix, markers)
    return prefix, markers, sample_ids


def check_grm(gcta, work_dir, algorithm, threads):
    prefix, markers, sample_ids = make_fixture(work_dir)
    output = run_gcta(
        gcta,
        work_dir,
        "grm",
        ["--bfile", prefix, "--make-grm-v1", "--make-grm-alg", algorithm],
        threads=threads,
    )
    expected_values, expected_counts = expected_legacy_grm(markers, algorithm)
    assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected_values, output)
    assert_floats_close(read_floats(output.with_suffix(".grm.N.bin")), expected_counts, output)
    expected_ids = "".join(f"{fid}\t{iid}\n" for fid, iid in sample_ids)
    assert output.with_suffix(".grm.id").read_text(encoding="utf-8") == expected_ids
    return output


def check_frequency(gcta, work_dir, threads):
    prefix, markers, _ = make_fixture(work_dir)
    output = run_gcta(gcta, work_dir, "freq", ["--bfile", prefix, "--freq-v1"], threads=threads)
    lines = [line.split() for line in output.with_suffix(".freq").read_text().splitlines()]
    assert len(lines) == len(markers)
    for index, (fields, marker) in enumerate(zip(lines, markers)):
        observed = [value for value in marker if value is not None]
        expected = sum(observed) / (2.0 * len(observed))
        assert fields[:2] == [f"rs{index}", "A"]
        assert len(fields) == 3
        assert math.isclose(float(fields[2]), expected, rel_tol=1e-12, abs_tol=1e-12)
    return output.with_suffix(".freq")


def main():
    parser = argparse.ArgumentParser(description="GCTA explicit legacy compatibility regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument(
        "--scenario",
        required=True,
        choices=("grm", "frequency", "threads_alias", "thread_conflict"),
    )
    args = parser.parse_args()
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)

    if args.scenario == "grm":
        for algorithm in (0, 1):
            single = check_grm(args.gcta, work / f"algorithm_{algorithm}" / "single", algorithm, 1)
            multi = check_grm(args.gcta, work / f"algorithm_{algorithm}" / "multi", algorithm, 4)
            assert_grm_equal(single, multi)
    elif args.scenario == "frequency":
        single = check_frequency(args.gcta, work / "single", 1)
        multi = check_frequency(args.gcta, work / "multi", 4)
        assert single.read_bytes() == multi.read_bytes()
    elif args.scenario == "threads_alias":
        prefix, _, _ = make_fixture(work / "input")
        standard = run_gcta(
            args.gcta, work, "standard", ["--bfile", prefix, "--make-grm"], threads=4
        )
        alias = run_gcta(
            args.gcta,
            work,
            "alias",
            ["--bfile", prefix, "--make-grm"],
            threads=4,
            thread_flag="--threads",
        )
        assert_grm_equal(standard, alias)
    else:
        prefix, _, _ = make_fixture(work / "input")
        expect_error(
            args.gcta,
            work,
            "thread_conflict",
            ["--bfile", prefix, "--make-grm", "--threads", 4, "--out", work / "result"],
        )
        console = (work / "thread_conflict.console").read_text(encoding="utf-8")
        assert "can't set the number of threads using both --thread-num and --threads" in console
    print(f"Passed explicit legacy regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
