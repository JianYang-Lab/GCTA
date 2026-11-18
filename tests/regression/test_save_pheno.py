#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink

SAMPLE_COUNT = 10
MARKERS = [
    [index % 3 for index in range(SAMPLE_COUNT)],
    [(2 * index) % 3 for index in range(SAMPLE_COUNT)],
    [0, 1, 2, 0, 1, 2, 0, 1, 2, 0],
    [2, 2, 1, 0, 0, 1, 2, 0, 1, 2],
]
COVARIATES = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
RESIDUALS = [0.2, -0.1, 0.3, 0.0, -0.4, 0.1, 0.2, -0.3, 0.4, -0.2]
MISSING_PHENOTYPE = 8


def make_fixture(work_dir):
    prefix = Path(work_dir) / "input" / "data"
    write_plink(prefix, MARKERS)
    phenotypes = [0.3 + 0.7 * covariate + residual for covariate, residual in zip(COVARIATES, RESIDUALS)]
    pheno = prefix.with_suffix(".phen")
    pheno.write_text(
        "".join(
            f"F1 I{index} {-9 if index == MISSING_PHENOTYPE else phenotypes[index]:.17g}\n"
            for index in range(SAMPLE_COUNT)
        ),
        encoding="utf-8",
    )
    qcovar = prefix.with_suffix(".qcovar")
    qcovar.write_text(
        "".join(f"F1 I{index} {COVARIATES[index]:.17g}\n" for index in range(SAMPLE_COUNT)),
        encoding="utf-8",
    )
    return prefix, pheno, qcovar, phenotypes


def expected_cphen(phenotypes):
    indices = [index for index in range(SAMPLE_COUNT) if index != MISSING_PHENOTYPE]
    covariates = [COVARIATES[index] for index in indices]
    observed = [phenotypes[index] for index in indices]
    mean_covariate = sum(covariates) / len(covariates)
    mean_phenotype = sum(observed) / len(observed)
    centered_covariates = [value - mean_covariate for value in covariates]
    slope = sum(
        covariate * (value - mean_phenotype)
        for covariate, value in zip(centered_covariates, observed)
    ) / sum(value * value for value in centered_covariates)
    intercept = mean_phenotype - slope * mean_covariate
    return [
        (index, phenotypes[index] - intercept - slope * COVARIATES[index])
        for index in indices
    ]


def check_output(gcta, work_dir, threads):
    prefix, pheno, qcovar, phenotypes = make_fixture(work_dir)
    output = run_gcta(
        gcta,
        work_dir,
        "fastgwa",
        [
            "--bfile", prefix,
            "--fastGWA-lr",
            "--pheno", pheno,
            "--qcovar", qcovar,
            "--save-pheno",
            "--seed", 12345,
        ],
        threads=threads,
    )
    rows = read_table(Path(f"{output}.fastGWA.cphen"))
    expected = expected_cphen(phenotypes)
    assert len(rows) == len(expected)
    for fields, (index, residual) in zip(rows, expected):
        assert fields[:2] == ["F1", f"I{index}"]
        assert math.isclose(float(fields[2]), residual, rel_tol=2e-6, abs_tol=2e-6)
    return output


def main():
    parser = argparse.ArgumentParser(description="GCTA adjusted phenotype output regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", choices=("save_pheno",), default="save_pheno")
    args = parser.parse_args()

    single = check_output(args.gcta, Path(args.work) / "single", 1)
    multi = check_output(args.gcta, Path(args.work) / "multi", 4)
    assert Path(f"{single}.fastGWA.cphen").read_bytes() == Path(f"{multi}.fastGWA.cphen").read_bytes()
    assert_tables_close(single.with_suffix(".fastGWA"), multi.with_suffix(".fastGWA"))
    print("Passed fastGWA adjusted phenotype output and thread regression checks")


if __name__ == "__main__":
    main()
