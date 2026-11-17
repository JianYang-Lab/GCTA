#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import read_table, run_gcta
from test_errors import expect_error


def beta_density(value, first, second):
    normalizer = math.gamma(first + second) / (math.gamma(first) * math.gamma(second))
    return normalizer * value ** (first - 1.0) * (1.0 - value) ** (second - 1.0)


def expected_acat(rows):
    weights = [
        (beta_density(af, 1.0, 25.0) / beta_density(af, 0.5, 0.5)) ** 2
        for af, _ in rows
    ]
    weight_sum = sum(weights)
    statistic = sum(
        weight / weight_sum * math.tan((0.5 - p_value) * math.pi)
        for weight, (_, p_value) in zip(weights, rows)
    )
    if abs(statistic) > 1.0:
        result = math.atan(1.0 / statistic) / math.pi
        return result if statistic > 0.0 else 1.0 + result
    return 0.5 - math.atan(statistic) / math.pi


def write_inputs(work):
    inputs = work / "input"
    inputs.mkdir(parents=True, exist_ok=True)
    genes = inputs / "genes.txt"
    snps = inputs / "snps.txt"
    genes.write_text("1 100 200 GENE1\n", encoding="utf-8")
    header = ["CHR", "SNP", "POS", "A1", "A2", "N", "AF", "BETA", "SE", "Z", "LOG10P", "EXTRA", "P"]
    included = [(0.02, 0.05), (0.03, 0.2)]
    rows = [
        ["1", "rs1", "150", "A", "G", "1000", "0.02", "0", "0", "0", "0", "0", "0.05"],
        ["1", "rs2", "160", "A", "G", "1000", "0.03", "0", "0", "0", "0", "0", "0.2"],
        ["1", "rs3", "170", "A", "G", "1000", "0.20", "0", "0", "0", "0", "0", "0.01"],
        ["1", "rs4", "500", "A", "G", "1000", "0.02", "0", "0", "0", "0", "0", "0.01"],
    ]
    snps.write_text(
        " ".join(header) + "\n" + "".join(" ".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    return genes, snps, expected_acat(included)


def main():
    parser = argparse.ArgumentParser(description="GCTA ACAT regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("output", "missing_input"))
    args = parser.parse_args()
    work = Path(args.work).resolve()
    genes, snps, expected = write_inputs(work)

    if args.scenario == "missing_input":
        expect_error(
            args.gcta,
            work,
            "missing_gene_list",
            ["--acat", "--gene-list", work / "absent.genes", "--snp-list", snps,
             "--out", work / "missing.csv"],
        )
    else:
        output = run_gcta(
            args.gcta,
            work,
            "acat.csv",
            ["--acat", "--gene-list", genes, "--snp-list", snps,
             "--max-maf", 0.05, "--min-mac", 20],
        )
        rows = read_table(output)
        assert rows[0] == ["CHR", "GENE", "START", "END", "SNP_NUM", "P_ACAT"]
        assert len(rows) == 2
        assert rows[1][:5] == ["1", "GENE1", "100", "200", "2"]
        assert math.isclose(float(rows[1][5]), expected, rel_tol=1e-6, abs_tol=1e-12), (
            f"ACAT P value differs: {rows[1][5]} != {expected}"
        )
    print(f"Passed ACAT regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
