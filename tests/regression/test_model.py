#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import assert_tables_close, run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes

SCENARIOS = ("model_lr", "model_mlm", "model_exact", "model_glmm")


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
	baseline = Path(args.data).parent / "baselines" / f"{args.scenario}.fastGWA"
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
