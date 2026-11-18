#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, run_gcta
from pgen import write_pgen
from plink import expected_grm, synthetic_markers
from test_geno import MARKERS, assert_frequencies, assert_genotypes


def check_grm(gcta, prefix, work_dir, markers, threads):
    output = run_gcta(gcta, work_dir, "grm", ["--pfile", prefix, "--make-grm"], threads=threads)
    expected_values, expected_counts = expected_grm(markers, 0)
    assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected_values, output)
    assert read_floats(output.with_suffix(".grm.N.bin")) == expected_counts
    return output


def check_hardcall(gcta, work_dir, threads):
    prefix = Path(work_dir) / "input" / "synthetic"
    write_pgen(prefix, MARKERS)
    grm = check_grm(gcta, prefix, work_dir, MARKERS, threads)
    freq = run_gcta(gcta, work_dir, "freq", ["--pfile", prefix, "--freq"], threads=threads)
    assert_frequencies(freq.with_suffix(".frq"))
    recode = run_gcta(gcta, work_dir, "recode", ["--pfile", prefix, "--recodet", "raw"], threads=threads)
    assert_genotypes(recode.with_suffix(".xmat"))
    return grm, freq.with_suffix(".frq"), recode.with_suffix(".xmat")


def main():
    parser = argparse.ArgumentParser(description="GCTA PGEN regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("hardcall", "dosage"))
    args = parser.parse_args()

    if args.scenario == "hardcall":
        single = check_hardcall(args.gcta, Path(args.work) / "single", 1)
        multi = check_hardcall(args.gcta, Path(args.work) / "multi", 4)
        assert_grm_equal(single[0], multi[0])
        assert single[1].read_bytes() == multi[1].read_bytes()
        assert single[2].read_bytes() == multi[2].read_bytes()
    else:
        for name, markers in (("complete", synthetic_markers()), ("missing", MARKERS)):
            prefix = Path(args.work) / name / "input" / "synthetic"
            write_pgen(prefix, markers, dosage=True)
            single = check_grm(args.gcta, prefix, Path(args.work) / name / "single", markers, 1)
            multi = check_grm(args.gcta, prefix, Path(args.work) / name / "multi", markers, 4)
            assert_grm_equal(single, multi)
    print(f"Passed PGEN regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
