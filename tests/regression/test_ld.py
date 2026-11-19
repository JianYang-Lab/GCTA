#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import struct

from gcta_runner import assert_floats_close, run_gcta
from plink import synthetic_markers, write_plink


def expected_values(markers, window, positions=None, chromosomes=None):
    standardized = []
    for marker in markers:
        observed = [value for value in marker if value is not None]
        frequency = sum(observed) / (2.0 * len(observed))
        variance = 2.0 * frequency * (1.0 - frequency)
        scale = math.sqrt(1.0 / variance) if variance >= 1e-50 else 0.0
        standardized.append([0.0 if value is None else (value - 2.0 * frequency) * scale for value in marker])

    expected = []
    sample_count = len(markers[0])
    for first in range(len(markers)):
        last = len(markers) if window else first + 1
        if positions is not None:
            last = first + 1
            while (last < len(markers) and chromosomes[last] == chromosomes[first]
                   and positions[last] - positions[first] <= window):
                last += 1
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


def check_output(gcta, prefix, work_dir, window, threads, expected, extra=()):
    output = run_gcta(
        gcta,
        work_dir,
        "ld",
        ["--bfile", prefix, "--ld-matrix", "--ld-wind", window, *extra],
        threads=threads,
    )
    actual = read_zld(output.with_suffix(".zld"))
    assert len(actual) == len(expected)
    for index, (actual_row, expected_row) in enumerate(zip(actual, expected)):
        assert_floats_close(actual_row, expected_row, f"{output}.zld row {index}")
    return output.with_suffix(".zld").read_bytes()


def check_boundaries(gcta, work):
    markers = synthetic_markers(sample_count=37, marker_count=289)
    markers[70][3] = None
    chromosomes = [1] * 281 + [2] * 8
    positions = [1000 + (i % 281) * 10 for i in range(len(markers))]
    prefix = work / "input" / "boundaries"
    write_plink(prefix, markers, chromosomes)
    bim = prefix.with_suffix(".bim")
    rows = [line.split() for line in bim.read_text().splitlines()]
    for row, position in zip(rows, positions):
        row[3] = str(position)
    bim.write_text("".join("\t".join(row) + "\n" for row in rows))
    excluded = {10, 101, 205, 280}
    exclusion = work / "exclude.txt"
    exclusion.write_text("".join(f"rs{i}\n" for i in sorted(excluded)))
    kept = [i for i in range(len(markers)) if i not in excluded]
    expected = expected_values([markers[i] for i in kept], 1000,
                               [positions[i] for i in kept], [chromosomes[i] for i in kept])
    outputs = [check_output(gcta, prefix, work / f"threads_{threads}", 1, threads, expected,
                            ("--exclude", exclusion)) for threads in (1, 4)]
    assert outputs[0] == outputs[1], "LD boundary fixture differs between thread counts"


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
        expected = expected_values(synthetic_markers(), window)
        single = check_output(args.gcta, prefix, work / f"window_{window}" / "single", window, 1, expected)
        multi = check_output(args.gcta, prefix, work / f"window_{window}" / "multi", window, 4, expected)
        assert single == multi, f"LD output differs between thread counts for window {window}"
    check_boundaries(args.gcta, work / "boundaries")
    print("Passed LD matrix numeric, window, and thread regression checks")


if __name__ == "__main__":
    main()
