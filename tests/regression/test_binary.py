#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import struct

from gcta_runner import read_table, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes

SCENARIOS = ("fastgwa_lr_savebin", "fastgwa_glmm_savebin")


def assert_close(actual, expected, label, abs_tol=1e-6):
    assert math.isclose(actual, expected, rel_tol=2e-5, abs_tol=abs_tol), (
        f"{label}: {actual} != {expected}"
    )


def run_scenario(gcta, work_dir, scenario, threads):
    markers, chromosomes = family_markers()
    prefix = Path(work_dir) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, qcovar = write_family_phenotypes(prefix, markers)
    option = "--fastGWA-lr" if scenario == "fastgwa_lr_savebin" else "--fastGWA-mlm-binary"
    arguments = [
        "--bfile", prefix,
        option,
        "--pheno", pheno,
        "--qcovar", qcovar,
        "--seed", 12345,
    ]
    if scenario == "fastgwa_glmm_savebin":
        sparse = run_gcta(
            gcta,
            work_dir,
            "sparse",
            ["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05],
            threads=threads,
        )
        arguments.extend(("--grm-sparse", sparse, "--mpheno", 2))
    text = run_gcta(gcta, work_dir, "text", arguments, threads=threads)
    binary = run_gcta(
        gcta, work_dir, "binary", [*arguments, "--save-bin"], threads=threads
    )
    text_path = Path(f"{text}.fastGWA")
    binary_path = Path(f"{binary}.fastGWA.bin")
    snpinfo_path = Path(f"{binary}.fastGWA.snpinfo")
    text_rows = read_table(text_path)
    snpinfo_rows = read_table(snpinfo_path)
    assert len(text_rows) == len(snpinfo_rows) == 2048
    assert snpinfo_rows[0] == ["CHR", "SNP", "POS", "A1", "A2"]
    for text_row, snpinfo_row in zip(text_rows[1:], snpinfo_rows[1:]):
        assert text_row[:5] == snpinfo_row

    layout = "<fffdI" if scenario == "fastgwa_lr_savebin" else "<fffddI"
    record = struct.Struct(layout)
    data = binary_path.read_bytes()
    assert len(data) == (len(text_rows) - 1) * record.size, (
        f"Truncated or trailing binary data: {binary_path}: {len(data)} bytes"
    )
    for index, values in enumerate(record.iter_unpack(data)):
        text_row = text_rows[index + 1]
        af, beta, se, p, count = values[:5]
        if scenario == "fastgwa_glmm_savebin":
            af, beta, se, p, padj, count = values
            assert_close(padj, float(text_row[12]), f"adjusted P at record {index}", 1e-12)
            assert_close(p, float(text_row[9]), f"unadjusted P at record {index}", 1e-12)
            beta_column, se_column = 10, 11
        else:
            assert_close(p, float(text_row[9]), f"P at record {index}", 1e-12)
            beta_column, se_column = 7, 8
        assert count == int(text_row[5])
        assert_close(af, float(text_row[6]), f"AF at record {index}")
        assert_close(beta, float(text_row[beta_column]), f"beta at record {index}")
        assert_close(se, float(text_row[se_column]), f"SE at record {index}")

    for output, count in ((text, 2047), (binary, 2047)):
        log = Path(f"{output}.log").read_text(encoding="utf-8")
        assert f"Saved {count} SNPs." in log
        assert "Error:" not in log
    return binary_path, snpinfo_path


def compare_baseline(binary, snpinfo, baseline, scenario):
    for actual, name in ((binary, "bin"), (snpinfo, "snpinfo")):
        expected = Path(baseline) / f"{scenario}.{name}"
        actual_data = actual.read_bytes()
        expected_data = expected.read_bytes()
        assert actual_data == expected_data, f"Output differs from baseline: {actual}"


def main():
    parser = argparse.ArgumentParser(description="GCTA binary output regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    parser.add_argument("--update-baselines", action="store_true")
    args = parser.parse_args()
    baseline = Path(args.data).parent / "baselines"
    if args.update_baselines:
        baseline.mkdir(parents=True, exist_ok=True)
    single = run_scenario(args.gcta, Path(args.work) / "single", args.scenario, 1)
    if args.update_baselines:
        for path in single:
            (baseline / f"{args.scenario}.{path.suffix[1:]}").write_bytes(path.read_bytes())
    multi = run_scenario(args.gcta, Path(args.work) / "multi", args.scenario, 4)
    compare_baseline(*single, baseline, args.scenario)
    compare_baseline(*multi, baseline, args.scenario)
    print(f"Passed binary output regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
