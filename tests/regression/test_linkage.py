#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import platform
import re
import subprocess

SYMBOLS = ("cblas_sgemm", "cblas_sgemv", "cblas_sdot", "dsyrk", "dgemm", "dpotrf", "LAPACKE_dsyev")


def run_with_bindings(gcta, work_dir, name, arguments):
	work_dir.mkdir(parents=True, exist_ok=True)
	environment = os.environ.copy()
	environment["LD_DEBUG"] = "bindings"
	environment["OMP_NUM_THREADS"] = "1"
	environment["MKL_NUM_THREADS"] = "1"
	output = work_dir / name
	command = [str(Path(gcta).resolve()), *map(str, arguments), "--out", output, "--thread-num", "1"]
	result = subprocess.run(
		command,
		cwd=work_dir,
		env=environment,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT,
		text=True,
		timeout=120,
		check=False,
	)
	(work_dir / f"{name}.bindings").write_text(result.stdout, encoding="utf-8")
	assert result.returncode == 0, f"Analysis failed during linkage audit: {name}"
	return result.stdout


def providers(output):
	result = {}
	pattern = re.compile(r"binding file .*? to (\S+) \[[0-9]+\]: normal symbol [`']([^' ]+)")
	for line in output.splitlines():
		match = pattern.search(line)
		if match and match.group(2) in SYMBOLS:
			result[match.group(2)] = match.group(1)
	return result


def main():
	parser = argparse.ArgumentParser(description="GCTA numerical linkage audit")
	parser.add_argument("--gcta", required=True)
	parser.add_argument("--data", required=True)
	parser.add_argument("--work", required=True)
	parser.add_argument("--scenario", default="numeric_symbols", choices=("numeric_symbols",))
	args = parser.parse_args()
	if platform.system() != "Linux":
		print("SKIP: this audit currently expects the Linux MKL backend")
		return
	gcta = Path(args.gcta).resolve()
	work = Path(args.work).resolve()
	legacy = Path(args.data).resolve() / "gcta_legacy" / "test"
	keep = work / "keep64.txt"
	work.mkdir(parents=True, exist_ok=True)
	with legacy.with_suffix(".fam").open(encoding="utf-8") as stream:
		keep.write_text("".join(stream.readlines()[:64]), encoding="utf-8")
	logs = []
	logs.append(run_with_bindings(gcta, work, "grm", ["--bfile", legacy, "--keep", keep, "--make-grm"]))
	grm = work / "grm"
	logs.append(run_with_bindings(gcta, work, "greml", [
		"--grm", grm, "--pheno", legacy.with_suffix(".phen"), "--reml",
		"--reml-priors", 0.5, "--reml-maxit", 100, "--reml-no-lrt",
	]))
	logs.append(run_with_bindings(gcta, work, "mlma", [
		"--bfile", legacy, "--keep", keep, "--mlma", "--pheno", legacy.with_suffix(".phen"),
		"--reml-priors", 0.5, "--reml-maxit", 100,
	]))
	logs.append(run_with_bindings(gcta, work, "pca", ["--grm", grm, "--pca", 10]))
	actual = providers("\n".join(logs))
	missing = [symbol for symbol in SYMBOLS if symbol not in actual]
	assert not missing, f"Symbols were not exercised by the linkage audit: {missing}"
	wrong = {symbol: library for symbol, library in actual.items() if "libmkl" not in library}
	assert not wrong, f"Numerical symbols resolve outside MKL: {wrong}"
	print("Passed numerical linkage audit: all required symbols use MKL")


if __name__ == "__main__":
	main()
