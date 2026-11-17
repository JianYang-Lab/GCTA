#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import assert_tables_close, read_table, run_gcta

SCENARIOS = ("greml", "pca", "fastgwa_lr", "mlma")


def write_keep_file(data_dir, work_dir):
	source = Path(data_dir) / "gcta_legacy" / "test.fam"
	keep = Path(work_dir) / "keep64.txt"
	keep.parent.mkdir(parents=True, exist_ok=True)
	with source.open(encoding="utf-8") as stream:
		ids = [next(stream).split()[:2] for _ in range(64)]
	keep.write_text("".join(f"{fid}\t{iid}\n" for fid, iid in ids), encoding="utf-8")
	return keep


def make_grm(gcta, prefix, keep, work_dir, threads):
	return run_gcta(
		gcta,
		work_dir,
		"grm",
		["--bfile", prefix, "--keep", keep, "--make-grm"],
		threads=threads,
	)


def run_workflow(gcta, data_dir, work_dir, scenario, threads):
	legacy = Path(data_dir) / "gcta_legacy" / "test"
	pheno = legacy.with_suffix(".phen")
	keep = write_keep_file(data_dir, work_dir)
	if scenario == "greml":
		grm = make_grm(gcta, legacy, keep, work_dir, threads)
		output = run_gcta(
			gcta,
			work_dir,
			"greml",
			[
				"--grm", grm,
				"--pheno", pheno,
				"--reml",
				"--reml-priors", 0.5,
				"--reml-maxit", 100,
				"--reml-no-lrt",
				"--reml-est-fix",
			],
			threads=threads,
		)
		assert "Log-likelihood ratio converged." in output.with_suffix(".log").read_text(encoding="utf-8")
		return {"greml.hsq": output.with_suffix(".hsq")}
	if scenario == "pca":
		grm = make_grm(gcta, legacy, keep, work_dir, threads)
		output = run_gcta(
			gcta, work_dir, "pca", ["--grm", grm, "--pca", 10], threads=threads
		)
		return {
			"pca.eigenval": output.with_suffix(".eigenval"),
			"pca.eigenvec": output.with_suffix(".eigenvec"),
		}
	if scenario == "fastgwa_lr":
		output = run_gcta(
			gcta,
			work_dir,
			"fastgwa",
			["--bfile", legacy, "--keep", keep, "--fastGWA-lr", "--pheno", pheno, "--seed", 12345],
			threads=threads,
		)
		log = output.with_suffix(".log").read_text(encoding="utf-8")
		assert "Performing fastGWA linear regression analysis..." in log
		return {"fastgwa_lr.fastGWA": output.with_suffix(".fastGWA")}

	output = run_gcta(
		gcta,
		work_dir,
		"mlma",
		[
			"--bfile", legacy,
			"--keep", keep,
			"--mlma",
			"--pheno", pheno,
			"--reml-priors", 0.5,
			"--reml-maxit", 100,
		],
		threads=threads,
	)
	assert "Log-likelihood ratio converged." in output.with_suffix(".log").read_text(encoding="utf-8")
	return {"mlma.mlma": output.with_suffix(".mlma")}


def normalized_eigenvectors(path):
	rows = read_table(path)
	if len(rows) < 2:
		return rows
	column_count = len(rows[0])
	for column in range(2, column_count):
		pivot = max(range(len(rows)), key=lambda row: abs(float(rows[row][column])))
		if float(rows[pivot][column]) < 0:
			for row in rows:
				row[column] = str(-float(row[column]))
	return rows


def assert_eigenvectors_close(actual, expected):
	actual_path = Path(actual)
	actual_normalized = normalized_eigenvectors(actual)
	expected_normalized = normalized_eigenvectors(expected)
	actual_temporary = actual_path.parent / f"{actual_path.name}.normalized"
	actual_temporary.write_text(
		"".join("\t".join(row) + "\n" for row in actual_normalized), encoding="utf-8"
	)
	expected_temporary = actual_path.parent / f"{actual_path.name}.expected.normalized"
	expected_temporary.write_text(
		"".join("\t".join(row) + "\n" for row in expected_normalized), encoding="utf-8"
	)
	assert_tables_close(actual_temporary, expected_temporary)


def compare_outputs(actual, baseline, scenario):
	for name, actual_path in actual.items():
		baseline_path = Path(baseline) / name
		assert baseline_path.exists(), f"Missing baseline: {baseline_path}"
		if scenario == "pca" and name.endswith(".eigenvec"):
			assert_eigenvectors_close(actual_path, baseline_path)
		else:
			assert_tables_close(actual_path, baseline_path)


def main():
	parser = argparse.ArgumentParser(description="GCTA official workflow regression tests")
	parser.add_argument("--gcta", required=True)
	parser.add_argument("--data", required=True)
	parser.add_argument("--work", required=True)
	parser.add_argument("--scenario", required=True, choices=SCENARIOS)
	parser.add_argument("--update-baselines", action="store_true")
	args = parser.parse_args()
	baseline = Path(args.data).parent / "baselines"
	if args.update_baselines:
		baseline.mkdir(parents=True, exist_ok=True)

	single = run_workflow(args.gcta, args.data, Path(args.work) / "single", args.scenario, 1)
	if args.update_baselines:
		for name, path in single.items():
			(baseline / name).write_bytes(path.read_bytes())
	multi = run_workflow(args.gcta, args.data, Path(args.work) / "multi", args.scenario, 4)
	compare_outputs(single, baseline, args.scenario)
	compare_outputs(multi, baseline, args.scenario)
	print(f"Passed workflow regression scenario: {args.scenario}")


if __name__ == "__main__":
	main()
