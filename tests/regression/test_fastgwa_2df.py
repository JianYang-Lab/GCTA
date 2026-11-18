#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import struct

from gcta_runner import read_table, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes

SCENARIOS = ("linear_sandwich", "linear_nosandwich", "mixed_sandwich", "mixed_nosandwich")


def write_envir(prefix, sample_count):
    path = Path(prefix).parent / "envir.txt"
    path.write_text(
        "".join(
            f"F1 I{index}\t{((index * 7) % 13 - 6) / 3.0:.17g}\n"
            for index in range(sample_count)
        ),
        encoding="utf-8",
    )
    return path


def make_arguments(gcta, work, scenario, markers, chromosomes):
    prefix = work / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, qcovar = write_family_phenotypes(prefix, markers)
    envir = write_envir(prefix, len(markers[0]))
    arguments = [
        "--bfile", prefix,
        "--pheno", pheno,
        "--qcovar", qcovar,
        "--envir", envir,
        "--seed", 12345,
    ]
    if scenario.startswith("mixed"):
        sparse = run_gcta(
            gcta,
            work,
            "sparse",
            ["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05],
        )
        arguments.extend(("--grm-sparse", sparse, "--fastGWA-mlm"))
    else:
        arguments.append("--fastGWA-lr")
    if scenario.endswith("nosandwich"):
        arguments.append("--noSandwich")
    return arguments


def check_pair(text, binary, label):
    rows = read_table(f"{text}.fastGWA")
    assert rows[0][0] == "CHR"
    assert "BETA_G_by_E" in rows[0] and "P_2df" in rows[0]
    snpinfo = read_table(f"{binary}.fastGWA.snpinfo")
    assert len(rows) == len(snpinfo)
    record = struct.Struct("<9f3dI")
    data = Path(f"{binary}.fastGWA.bin").read_bytes()
    assert len(data) == (len(rows) - 1) * record.size, f"Unexpected binary size: {binary}"
    for (text_row, snp_row), values in zip(zip(rows[1:], snpinfo[1:]), record.iter_unpack(data)):
        assert text_row[:5] == snp_row
        af, beta_geno, beta_interaction, se_geno, se_interaction, cov, score_geno, \
            score_interaction, score, p_geno, p_interaction, p, count = values
        assert count == int(text_row[5]), f"{label}: N mismatch"
        numeric = (
            (af, text_row[6]), (beta_geno, text_row[7]), (beta_interaction, text_row[8]),
            (se_geno, text_row[9]), (se_interaction, text_row[10]), (cov, text_row[11]),
            (score_geno, text_row[12]), (score_interaction, text_row[13]), (score, text_row[14]),
            (p_geno, text_row[15]), (p_interaction, text_row[16]), (p, text_row[17]),
        )
        for observed, reference in numeric:
            reference = float(reference)
            if math.isnan(observed) and math.isnan(reference):
                continue
            assert math.isclose(observed, reference, rel_tol=2e-5, abs_tol=1e-8), (
                f"{label}: {observed} != {reference}"
            )
    return Path(f"{binary}.fastGWA.bin").read_bytes()


def run_scenario(gcta, work, scenario, threads):
    markers, chromosomes = family_markers()
    arguments = make_arguments(gcta, work, scenario, markers, chromosomes)
    text = run_gcta(gcta, work, "text", arguments, threads=threads)
    binary = run_gcta(gcta, work, "binary", [*arguments, "--save-bin"], threads=threads)
    return check_pair(text, binary, f"{scenario}/threads={threads}")


def main():
    parser = argparse.ArgumentParser(description="GCTA batched fastGWA 2df regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    args = parser.parse_args()
    work = Path(args.work)
    single = run_scenario(args.gcta, work / "single", args.scenario, 1)
    multi = run_scenario(args.gcta, work / "multi", args.scenario, 4)
    assert single == multi, f"Save-bin output differs between thread counts: {args.scenario}"
    print(f"Passed batched fastGWA 2df regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
