#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import assert_tables_close, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes

SCENARIOS = ("model_lr", "model_lr_reordered", "model_mlm", "model_exact", "model_glmm")


def run_roundtrip(gcta, work_dir, scenario, threads):
    markers, chromosomes = family_markers()
    prefix = Path(work_dir) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, qcovar = write_family_phenotypes(prefix, markers)
    arguments = [
        "--bfile", prefix,
        "--pheno", pheno,
        "--qcovar", qcovar,
        "--joint-covar",
        "--seed", 12345,
    ]
    if scenario in ("model_mlm", "model_exact", "model_glmm"):
        sparse = run_gcta(
            gcta,
            work_dir,
            "sparse",
            ["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05],
            threads=threads,
        )
        arguments.extend(("--grm-sparse", sparse))
    option = {
        "model_lr": "--fastGWA-lr",
        "model_lr_reordered": "--fastGWA-lr",
        "model_mlm": "--fastGWA-mlm",
        "model_exact": "--fastGWA-mlm-exact",
        "model_glmm": "--fastGWA-mlm-binary",
    }[scenario]
    arguments.append(option)
    if scenario == "model_glmm":
        arguments.extend(("--mpheno", 2))
    direct = run_gcta(gcta, work_dir, "direct", arguments, threads=threads)
    model = run_gcta(
        gcta, work_dir, "model", [*arguments, "--model-only"], threads=threads
    )
    load_prefix = prefix
    if scenario == "model_lr_reordered":
        load_prefix = Path(work_dir) / "input" / "family_reversed"
        write_plink(
            load_prefix,
            [list(reversed(marker)) for marker in markers],
            chromosomes,
        )
        fam_lines = prefix.with_suffix(".fam").read_text(encoding="utf-8").splitlines()
        load_prefix.with_suffix(".fam").write_text(
            "".join(f"{line}\n" for line in reversed(fam_lines)),
            encoding="utf-8",
        )
        direct_reordered = run_gcta(
            gcta,
            work_dir,
            "direct_reordered",
            ["--bfile", load_prefix, *arguments[2:]],
            threads=threads,
        )
        assert_tables_close(Path(f"{direct}.fastGWA"), Path(f"{direct_reordered}.fastGWA"))
    loaded = run_gcta(
        gcta,
        work_dir,
        "loaded",
        ["--bfile", load_prefix, "--load-model", Path(f"{model}.fastGWA")],
        threads=threads,
    )
    return Path(f"{direct}.fastGWA"), Path(f"{loaded}.fastGWA")


def main():
    parser = argparse.ArgumentParser(description="GCTA model roundtrip regression test")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    parser.add_argument("--update-baselines", action="store_true")
    args = parser.parse_args()
    baseline_scenario = "model_lr" if args.scenario == "model_lr_reordered" else args.scenario
    baseline = Path(args.data).parent / "baselines" / f"{baseline_scenario}.fastGWA"
    single = run_roundtrip(args.gcta, Path(args.work) / "single", args.scenario, 1)
    if args.update_baselines:
        baseline.write_bytes(single[0].read_bytes())
    multi = run_roundtrip(args.gcta, Path(args.work) / "multi", args.scenario, 4)
    for direct, loaded in (single, multi):
        assert_tables_close(direct, loaded)
        assert_tables_close(direct, baseline)
        assert_tables_close(loaded, baseline)
    print(f"Passed model roundtrip regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
