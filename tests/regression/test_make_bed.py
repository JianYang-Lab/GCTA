#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from bgen import write_bgen
from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, read_table, run_gcta
from plink import expected_grm, write_plink

MARKERS = [
    [2, None, 1, 0, 2, 0, 1, 2],
    [None, 1, 0, 2, None, 1, 2, 0],
    [1, 0, 2, None, 1, 2, 0, 1],
    [0, 2, None, 1, 0, 2, 1, 0],
    [0, 0, 1, 1, 2, 2, 1, 0],
]
KEEP_SAMPLES = (0, 2, 4, 5, 7)
KEEP_MARKERS = (0, 2, 4)


def make_fixture(work_dir, input_format):
    prefix = Path(work_dir) / "input" / "source"
    if input_format == "plink":
        write_plink(prefix, MARKERS)
        arguments = ["--bfile", prefix]
    else:
        write_bgen(prefix, MARKERS)
        arguments = ["--bgen", prefix.with_suffix(".bgen"), "--sample", prefix.with_suffix(".sample")]

    keep = Path(work_dir) / "input" / "keep.txt"
    keep.write_text("".join(f"F1 I{index}\n" for index in KEEP_SAMPLES), encoding="utf-8")
    exclude = Path(work_dir) / "input" / "exclude.txt"
    exclude.write_text("rs1\nrs3\n", encoding="utf-8")
    arguments.extend(("--keep", keep, "--exclude", exclude))

    expected_markers = [
        [MARKERS[index][sample] for sample in KEEP_SAMPLES]
        for index in KEEP_MARKERS
    ]
    expected_bed = Path(work_dir) / "input" / "expected"
    write_plink(expected_bed, expected_markers)
    return arguments, expected_bed.with_suffix(".bed"), expected_markers


def check_conversion(gcta, work_dir, input_format, threads):
    arguments, expected_bed, expected_markers = make_fixture(work_dir, input_format)
    converted = run_gcta(
        gcta,
        work_dir,
        "converted",
        [*arguments, "--make-bed"],
        threads=threads,
    )
    assert converted.with_suffix(".bed").read_bytes() == expected_bed.read_bytes()

    bim = read_table(converted.with_suffix(".bim"))
    assert bim == [["1", f"rs{index}", "0", str(1000 + index), "A", "G"] for index in KEEP_MARKERS]
    fam = read_table(converted.with_suffix(".fam"))
    assert [fields[:5] for fields in fam] == [["F1", f"I{index}", "0", "0", "0"] for index in KEEP_SAMPLES]
    assert all(math.isnan(float(fields[5])) for fields in fam)

    grm = run_gcta(
        gcta,
        work_dir,
        "roundtrip",
        ["--bfile", converted, "--make-grm"],
        threads=threads,
    )
    expected_values, expected_counts = expected_grm(expected_markers, 0)
    assert_floats_close(read_floats(grm.with_suffix(".grm.bin")), expected_values, grm)
    assert read_floats(grm.with_suffix(".grm.N.bin")) == expected_counts
    expected_ids = "".join(f"F1\tI{index}\n" for index in KEEP_SAMPLES)
    assert grm.with_suffix(".grm.id").read_text(encoding="utf-8") == expected_ids
    return converted, grm


def main():
    parser = argparse.ArgumentParser(description="GCTA PLINK binary output regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("plink", "bgen"))
    args = parser.parse_args()

    single = check_conversion(args.gcta, Path(args.work) / "single", args.scenario, 1)
    multi = check_conversion(args.gcta, Path(args.work) / "multi", args.scenario, 4)
    for suffix in ("bed", "bim", "fam"):
        assert single[0].with_suffix(f".{suffix}").read_bytes() == multi[0].with_suffix(f".{suffix}").read_bytes()
    assert_grm_equal(single[1], multi[1])
    print(f"Passed PLINK binary conversion scenario: {args.scenario}")


if __name__ == "__main__":
    main()
