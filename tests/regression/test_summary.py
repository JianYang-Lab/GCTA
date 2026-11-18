#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes

MARKER_INDICES = (1, 2)


def normal_pvalue(z_value):
    return math.erfc(abs(z_value) / math.sqrt(2.0))


def make_fixture(work_dir):
    markers, chromosomes = family_markers()
    prefix = Path(work_dir) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, _ = write_family_phenotypes(prefix, markers)
    phenotypes = [float(line.split()[2]) for line in pheno.read_text().splitlines()]
    statistics = []
    for index in MARKER_INDICES:
        genotypes = [float(value) for value in markers[index]]
        marker_mean = sum(genotypes) / len(genotypes)
        phenotype_mean = sum(phenotypes) / len(phenotypes)
        covariance = sum(
            (genotype - marker_mean) * (phenotype - phenotype_mean)
            for genotype, phenotype in zip(genotypes, phenotypes)
        )
        variance = sum((genotype - marker_mean) ** 2 for genotype in genotypes)
        beta = covariance / variance
        residuals = [
            phenotype - phenotype_mean - beta * (genotype - marker_mean)
            for genotype, phenotype in zip(genotypes, phenotypes)
        ]
        se = math.sqrt(sum(value * value for value in residuals) / (len(genotypes) - 2) / variance)
        z_value = beta / se
        statistics.append({
            "index": index,
            "beta": beta,
            "se": se,
            "af": marker_mean / 2.0,
            "p": normal_pvalue(z_value),
            "chisq": z_value * z_value,
        })
    meta = Path(work_dir) / "meta.ma"
    meta.write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{stat['index']} A G {stat['af']} {stat['beta']} {stat['se']} {stat['p']} {len(phenotypes)}\n"
            for stat in statistics
        ),
        encoding="utf-8",
    )
    genes = Path(work_dir) / "genes.txt"
    genes.write_text("1 1001 1002 GENE1\n", encoding="utf-8")
    return prefix, meta, genes, statistics


def run_scenario(gcta, work_dir, scenario, threads):
    prefix, meta, genes, statistics = make_fixture(work_dir)
    if scenario == "cojo_joint":
        snplist = Path(work_dir) / "rs.snplist"
        snplist.write_text("rs1\n", encoding="utf-8")
        output = run_gcta(
            gcta,
            work_dir,
            "joint",
            ["--bfile", prefix, "--extract", snplist, "--cojo-file", meta, "--cojo-joint"],
            threads=threads,
        )
        rows = read_table(Path(f"{output}.jma.cojo")
        )
        assert len(rows) == 2
        expected = statistics[0]
        assert rows[1][1] == "rs1"
        assert math.isclose(float(rows[1][5]), expected["beta"], rel_tol=1e-5, abs_tol=1e-6)
        assert math.isclose(float(rows[1][6]), expected["se"], rel_tol=1e-5, abs_tol=1e-6)
        assert math.isclose(float(rows[1][10]), expected["beta"], rel_tol=1e-5, abs_tol=1e-6)
        assert math.isclose(float(rows[1][11]), expected["se"], rel_tol=2e-5, abs_tol=2e-6)
        assert math.isclose(float(rows[1][12]), expected["p"], rel_tol=2e-4, abs_tol=1e-6)
        return Path(f"{output}.jma.cojo")

    option = "--fastBAT" if scenario == "fastbat" else "--mBAT-combo"
    output = run_gcta(
        gcta,
        work_dir,
        scenario,
        ["--bfile", prefix, option, meta, "--fastBAT-gene-list" if scenario == "fastbat" else "--mBAT-gene-list", genes],
        threads=threads,
    )
    suffix = ".gene.fastbat" if scenario == "fastbat" else ".gene.assoc.mbat"
    result = Path(f"{output}{suffix}")
    rows = read_table(result)
    assert len(rows) == 2 and rows[1][0] == "GENE1"
    best = min(statistics, key=lambda stat: stat["p"])
    assert rows[1][4] == "2"
    if scenario == "fastbat":
        assert math.isclose(float(rows[1][7]), sum(stat["chisq"] for stat in statistics), rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(float(rows[1][9]), best["p"], rel_tol=2e-5, abs_tol=1e-6)
        assert rows[1][10] == f"rs{best['index']}"
        assert 0.0 < float(rows[1][8]) < 1.0
    else:
        assert rows[1][7] == f"rs{best['index']}"
        assert math.isclose(float(rows[1][8]), best["p"], rel_tol=2e-5, abs_tol=1e-6)
        assert rows[1][9] == "2"
        assert float(rows[1][10]) > 0.0
        assert 0.0 < float(rows[1][11]) < 1.0
    return result


def main():
    parser = argparse.ArgumentParser(description="GCTA summary-statistics workflow regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("cojo_joint", "fastbat", "mbat"))
    args = parser.parse_args()
    single = run_scenario(args.gcta, Path(args.work) / "single", args.scenario, 1)
    multi = run_scenario(args.gcta, Path(args.work) / "multi", args.scenario, 4)
    assert_tables_close(single, multi)
    print(f"Passed summary-statistics workflow scenario: {args.scenario}")


if __name__ == "__main__":
    main()
