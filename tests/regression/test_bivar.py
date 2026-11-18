#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes
from test_errors import expect_error


def make_fixture(work_dir, gcta):
    markers, chromosomes = family_markers()
    prefix = Path(work_dir) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, _ = write_family_phenotypes(prefix, markers)
    rows = []
    for index, line in enumerate(pheno.read_text(encoding="utf-8").splitlines()):
        fields = line.split()
        first = float(fields[2])
        second = 0.4 * first + 0.1 * (index % 7 - 3)
        fields[2] = str(first)
        fields.append(str(second))
        rows.append("\t".join(fields) + "\n")
    pheno_path = Path(work_dir) / "two.phen"
    pheno_path.write_text("".join(rows), encoding="utf-8")
    grm = run_gcta(gcta, work_dir, "grm", ["--bfile", prefix, "--make-grm"], threads=1)
    return grm, pheno_path


def run_bivar(gcta, work_dir, threads):
    grm, pheno = make_fixture(work_dir, gcta)
    output = run_gcta(
        gcta,
        work_dir,
        "bivar",
        ["--grm", grm, "--pheno", pheno, "--reml-bivar", "--reml-maxit", 100, "--reml-no-lrt"],
        threads=threads,
    )
    log = output.with_suffix(".log").read_text(encoding="utf-8")
    assert "Log-likelihood ratio converged." in log
    values = {}
    for row in read_table(output.with_suffix(".hsq")):
        if len(row) < 2:
            continue
        try:
            values[row[0]] = float(row[1])
        except ValueError:
            continue
    for trait in (1, 2):
        assert math.isclose(
            values[f"V(G)_tr{trait}"] + values[f"V(e)_tr{trait}"],
            values[f"Vp_tr{trait}"],
            rel_tol=1e-5,
            abs_tol=1e-6,
        )
        assert math.isclose(
            values[f"V(G)_tr{trait}"] / values[f"Vp_tr{trait}"],
            values[f"V(G)/Vp_tr{trait}"],
            rel_tol=1e-5,
            abs_tol=1e-6,
        )
    assert -1.0 <= values["rG"] <= 1.0
    assert values["n"] == 256
    assert math.isfinite(values["logL"])
    return output.with_suffix(".hsq")


def main():
    parser = argparse.ArgumentParser(description="GCTA bivariate REML regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("reml", "bad_priors"))
    args = parser.parse_args()
    work = Path(args.work)

    if args.scenario == "reml":
        single = run_bivar(args.gcta, work / "single", 1)
        multi = run_bivar(args.gcta, work / "multi", 4)
        assert_tables_close(single, multi)
    else:
        grm, pheno = make_fixture(work, args.gcta)
        expect_error(
            args.gcta,
            work,
            "bad_priors",
            ["--grm", grm, "--pheno", pheno, "--reml-bivar",
             "--reml-priors", 0.5, 0.5, 0.5, "--reml-no-lrt", "--out", work / "result"],
        )
        console = (work / "bad_priors.console").read_text(encoding="utf-8")
        assert "At least 5 prior values should be specified" in console
    print(f"Passed bivariate REML scenario: {args.scenario}")


if __name__ == "__main__":
    main()
