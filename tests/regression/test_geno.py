#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import run_gcta
from plink import write_plink

MARKERS = [
    [2, None, 1, 0, 2],
    [None, 1, 0, 2, None],
    [1, 0, 2, None, 1],
    [0, 2, None, 1, 0],
    [0, 0, 0, 0, 0],
    [2, 2, 2, 2, 2],
]


def marker_statistics(marker):
    observed = [value for value in marker if value is not None]
    frequency = sum(observed) / (2.0 * len(observed))
    return frequency, 2 * len(observed), 2.0 * frequency * (1.0 - frequency)


def assert_metadata(fields, marker_index, frequency, chromosome_count):
    assert fields[:5] == ["1", f"rs{marker_index}", str(1000 + marker_index), "A", "G"]
    assert math.isclose(float(fields[5]), frequency, rel_tol=0.0, abs_tol=1e-6)
    assert fields[6] == str(chromosome_count)


def assert_frequencies(path):
    lines = [line.split("\t") for line in Path(path).read_text(encoding="utf-8").splitlines()]
    assert lines[0] == ["CHR", "SNP", "POS", "A1", "A2", "AF", "NCHROBS"]
    assert len(lines) == len(MARKERS) + 1
    for marker_index, marker in enumerate(MARKERS):
        frequency, chromosome_count, _ = marker_statistics(marker)
        assert_metadata(lines[marker_index + 1], marker_index, frequency, chromosome_count)


def assert_genotypes(path):
    lines = [line.split("\t") for line in Path(path).read_text(encoding="utf-8").splitlines()]
    sample_count = len(MARKERS[0])
    assert lines[0][:7] == ["CHR", "SNP", "POS", "A1", "A2", "AF", "NCHROBS"]
    assert lines[0][7:] == [f"F1|I{index}" for index in range(sample_count)]
    expected_markers = [
        (index, marker) for index, marker in enumerate(MARKERS) if marker_statistics(marker)[2] >= 1e-50
    ]
    assert len(lines) == len(expected_markers) + 1
    for fields, (marker_index, marker) in zip(lines[1:], expected_markers):
        frequency, chromosome_count, _ = marker_statistics(marker)
        assert_metadata(fields, marker_index, frequency, chromosome_count)
        expected = ["NA" if value is None else str(value) for value in marker]
        assert fields[7:] == expected, f"Genotype mismatch for rs{marker_index}"


def main():
    parser = argparse.ArgumentParser(description="GCTA genotype decoding regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", choices=("decoding",), default="decoding")
    args = parser.parse_args()
    prefix = Path(args.work) / "input" / "decoding"
    write_plink(prefix, MARKERS)

    for threads in (1, 4):
        output = run_gcta(
            args.gcta,
            Path(args.work) / f"threads_{threads}",
            "freq",
            ["--bfile", prefix, "--freq"],
            threads=threads,
        )
        assert_frequencies(output.with_suffix(".frq"))
        recode = run_gcta(
            args.gcta,
            Path(args.work) / f"threads_{threads}",
            "recode",
            ["--bfile", prefix, "--recodet", "raw"],
            threads=threads,
        )
        assert_genotypes(recode.with_suffix(".xmat"))

    first = Path(args.work) / "threads_1"
    fourth = Path(args.work) / "threads_4"
    assert (first / "freq.frq").read_bytes() == (fourth / "freq.frq").read_bytes()
    assert (first / "recode.xmat").read_bytes() == (fourth / "recode.xmat").read_bytes()
    print("Passed genotype decoding, frequency, missing-value, and thread regression checks")


if __name__ == "__main__":
    main()
