#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import struct

from gcta_runner import assert_floats_close, run_gcta
from plink import synthetic_markers, write_plink


def expected_values(markers, window):
    standardized = []
    for marker in markers:
        frequency = sum(marker) / (2.0 * len(marker))
        variance = 2.0 * frequency * (1.0 - frequency)
        scale = math.sqrt(1.0 / variance) if variance >= 1e-50 else 0.0
        standardized.append([(value - 2.0 * frequency) * scale for value in marker])

    expected = []
    sample_count = len(markers[0])
    for first in range(len(markers)):
        last = len(markers) if window else first + 1
        expected.append([
            sum(standardized[first][sample] * standardized[second][sample]
                for sample in range(sample_count)) / (sample_count - 1)
            for second in range(first, last)
        ])
    return expected


def read_zld(path):
    data = Path(path).read_bytes()
    offset = 0
    rows = []
    while offset < len(data):
        count = struct.unpack_from("<I", data, offset)[0]
        offset += 4
        rows.append(list(struct.unpack_from(f"<{count}f", data, offset)))
        offset += 4 * count
    assert offset == len(data), f"Truncated or trailing LD data: {path}"
    return rows


def check_output(gcta, prefix, work_dir, window, threads):
    output = run_gcta(
        gcta,
        work_dir,
        "ld",
        ["--bfile", prefix, "--ld-matrix", "--ld-wind", window],
        threads=threads,
    )
    actual = read_zld(output.with_suffix(".zld"))
    expected = expected_values(synthetic_markers(), window)
    assert len(actual) == len(expected)
    for index, (actual_row, expected_row) in enumerate(zip(actual, expected)):
        assert_floats_close(actual_row, expected_row, f"{output}.zld row {index}")
    return output.with_suffix(".zld").read_bytes()


def main():
    parser = argparse.ArgumentParser(description="GCTA LD matrix regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", default="matrix", choices=("matrix",))
    args = parser.parse_args()
    work = Path(args.work)
    prefix = work / "input" / "synthetic"
    write_plink(prefix, synthetic_markers())

    for window in (0, 1000):
        single = check_output(args.gcta, prefix, work / f"window_{window}" / "single", window, 1)
        multi = check_output(args.gcta, prefix, work / f"window_{window}" / "multi", window, 4)
        assert single == multi, f"LD output differs between thread counts for window {window}"
    print("Passed LD matrix numeric, window, and thread regression checks")


if __name__ == "__main__":
    main()
