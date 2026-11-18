#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from bgen import write_bgen
from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, run_gcta
from pgen import write_pgen
from plink import expected_grm, synthetic_markers, write_plink

MARKERS = synthetic_markers(sample_count=8, marker_count=8)
PARTITIONS = (MARKERS[:3], MARKERS[3:])


def make_fixture(work_dir, scenario):
    work_dir = Path(work_dir)
    prefixes = []
    for index, markers in enumerate(PARTITIONS):
        prefix = work_dir / "input" / f"part{index}"
        if scenario == "mbfile":
            write_plink(prefix, markers)
        elif scenario == "mbgen":
            write_bgen(prefix, markers)
        else:
            write_pgen(prefix, markers)
        prefixes.append(prefix)
    list_path = work_dir / "files.txt"
    if scenario == "mbgen":
        list_path.write_text("".join(f"{prefix.with_suffix('.bgen')}\n" for prefix in prefixes), encoding="utf-8")
    else:
        list_path.write_text("".join(f"{prefix}\n" for prefix in prefixes), encoding="utf-8")
    return prefixes, list_path


def main():
    parser = argparse.ArgumentParser(description="GCTA multi-file input regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("mbfile", "mbgen", "mpfile"))
    args = parser.parse_args()
    if args.scenario == "mbgen":
        valid_marker_count = sum(min(marker) != max(marker) for marker in MARKERS)
        expected_values = []
        for row in range(len(MARKERS[0])):
            for column in range(row + 1):
                total = 0.0
                for marker in MARKERS:
                    mean = sum(marker) / len(marker)
                    variance = sum((value - mean) ** 2 for value in marker) / (len(marker) - 1)
                    scale = math.sqrt(1.0 / variance) if variance > 1e-50 else 0.0
                    total += (marker[row] - mean) * (marker[column] - mean) * scale * scale
                expected_values.append(total / valid_marker_count)
        expected_counts = [float(valid_marker_count)] * 36
    else:
        expected_values, expected_counts = expected_grm(MARKERS, 0)
    outputs = []
    for threads in (1, 4):
        work = Path(args.work) / f"threads_{threads}"
        prefixes, list_path = make_fixture(work, args.scenario)
        arguments = [f"--{args.scenario}", list_path, "--make-grm"]
        if args.scenario == "mbgen":
            arguments.extend(("--sample", prefixes[0].with_suffix(".sample")))
        output = run_gcta(args.gcta, work, "grm", arguments, threads=threads)
        assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected_values, output)
        assert read_floats(output.with_suffix(".grm.N.bin")) == expected_counts
        outputs.append(output)
    assert_grm_equal(*outputs)
    print(f"Passed multi-file regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
