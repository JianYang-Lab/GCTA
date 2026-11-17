#!/usr/bin/env python3

import argparse
from pathlib import Path

from bgen import write_bgen
from gcta_runner import assert_grm_equal, run_gcta
from plink import synthetic_markers


def make_fixture(work_dir):
    markers = synthetic_markers(sample_count=8, marker_count=17)
    markers[2][7] = None
    markers[3][2] = None
    markers[4][1] = None
    markers[5][4] = None
    keep_indices = (1, 2, 4, 7)
    compact_markers = [[marker[index] for index in keep_indices] for marker in markers]
    full_prefix = Path(work_dir) / "input" / "full"
    write_bgen(full_prefix, markers)
    compact_prefix = Path(work_dir) / "input" / "compact"
    write_bgen(compact_prefix, compact_markers, sample_indices=keep_indices)
    keep = Path(work_dir) / "input" / "keep.txt"
    keep.write_text(
        "".join(f"F1\tI{index}\n" for index in keep_indices), encoding="utf-8"
    )
    return full_prefix, compact_prefix, keep


def run_pair(gcta, work_dir, threads):
    full_prefix, compact_prefix, keep = make_fixture(work_dir)
    subset = run_gcta(
        gcta,
        work_dir,
        "subset",
        ["--bgen", full_prefix.with_suffix(".bgen"), "--sample", full_prefix.with_suffix(".sample"), "--keep", keep, "--make-grm"],
        threads=threads,
    )
    compact = run_gcta(
        gcta,
        work_dir,
        "compact",
        ["--bgen", compact_prefix.with_suffix(".bgen"), "--sample", compact_prefix.with_suffix(".sample"), "--make-grm"],
        threads=threads,
    )
    assert_grm_equal(subset, compact)
    return subset


def main():
    parser = argparse.ArgumentParser(description="GCTA genotype format regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", default="bgen_subset", choices=("bgen_subset",))
    parser.add_argument("--update-baselines", action="store_true")
    args = parser.parse_args()
    baseline = Path(args.data).parent / "baselines" / "bgen_subset"
    single = run_pair(args.gcta, Path(args.work) / "single", 1)
    if args.update_baselines:
        for suffix in ("grm.bin", "grm.N.bin", "grm.id"):
            Path(f"{baseline}.{suffix}").write_bytes(Path(f"{single}.{suffix}").read_bytes())
    multi = run_pair(args.gcta, Path(args.work) / "multi", 4)
    for output in (single, multi):
        assert_grm_equal(output, baseline)
    print("Passed genotype format regression scenario: bgen_subset")


if __name__ == "__main__":
    main()
