#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import read_floats, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes


def solve_ols(grm_prefix, pheno_path):
    ids = [line.split() for line in Path(f"{grm_prefix}.grm.id").read_text().splitlines()]
    phenotypes = [float(line.split()[2]) for line in Path(pheno_path).read_text().splitlines()]
    assert len(ids) == len(phenotypes)
    mean = sum(phenotypes) / len(phenotypes)
    scale = math.sqrt(sum((value - mean) ** 2 for value in phenotypes) / (len(phenotypes) - 1))
    values = [(value - mean) / scale for value in phenotypes]
    grm = read_floats(Path(f"{grm_prefix}.grm.bin"))
    lhs00 = lhs01 = lhs11 = 0.0
    rhs_cp = [0.0, 0.0]
    rhs_sd = [0.0, 0.0]
    index = 0
    for row in range(len(values)):
        for column in range(row):
            relatedness = grm[index]
            index += 1
            lhs00 += 1.0
            lhs01 += relatedness
            lhs11 += relatedness * relatedness
            cross_product = values[row] * values[column]
            squared_difference = (values[row] - values[column]) ** 2
            rhs_cp[0] += cross_product
            rhs_cp[1] += relatedness * cross_product
            rhs_sd[0] += squared_difference
            rhs_sd[1] += relatedness * squared_difference
        index += 1
    lhs = [[lhs00, lhs01], [lhs01, lhs11]]
    determinant = lhs[0][0] * lhs[1][1] - lhs[0][1] ** 2
    inverse = [
        [lhs[1][1] / determinant, -lhs[0][1] / determinant],
        [-lhs[0][1] / determinant, lhs[0][0] / determinant],
    ]
    beta_cp = [
        inverse[0][0] * rhs_cp[0] + inverse[0][1] * rhs_cp[1],
        inverse[1][0] * rhs_cp[0] + inverse[1][1] * rhs_cp[1],
    ]
    beta_sd = [
        inverse[0][0] * rhs_sd[0] + inverse[0][1] * rhs_sd[1],
        inverse[1][0] * rhs_sd[0] + inverse[1][1] * rhs_sd[1],
    ]
    return beta_cp, beta_sd


def parse_output(path):
    sections = {}
    current = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line in ("HE-CP", "HE-SD"):
            current = line
            sections[current] = {}
        elif current and line.startswith("Intercept"):
            sections[current]["Intercept"] = float(line.split()[1])
        elif current and line.startswith("V(G)/Vp"):
            sections[current]["V(G)/Vp"] = float(line.split()[1])
    return sections


def main():
    parser = argparse.ArgumentParser(description="GCTA Haseman-Elston regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", default="output", choices=("output",))
    args = parser.parse_args()
    markers, chromosomes = family_markers()
    prefix = Path(args.work) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, _ = write_family_phenotypes(prefix, markers)
    grm = run_gcta(args.gcta, args.work, "grm", ["--bfile", prefix, "--make-grm"])
    expected_cp, _ = solve_ols(grm, pheno)

    outputs = []
    for threads in (1, 4):
        output = run_gcta(
            args.gcta,
            Path(args.work) / f"threads_{threads}",
            "he",
            ["--grm", grm, "--pheno", pheno, "--HEreg"],
            threads=threads,
        )
        sections = parse_output(output.with_suffix(".HEreg"))
        assert math.isclose(sections["HE-CP"]["V(G)/Vp"], expected_cp[1], rel_tol=2e-5, abs_tol=1e-6)
        assert math.isclose(sections["HE-CP"]["Intercept"], expected_cp[0], rel_tol=2e-5, abs_tol=1e-6)
        assert math.isfinite(sections["HE-SD"]["V(G)/Vp"])
        assert math.isfinite(sections["HE-SD"]["Intercept"])
        outputs.append(output.with_suffix(".HEreg"))
    assert outputs[0].read_bytes() == outputs[1].read_bytes()
    print("Passed HE-CP independent OLS oracle and HE-SD structure with 1/4-thread agreement")


if __name__ == "__main__":
    main()
