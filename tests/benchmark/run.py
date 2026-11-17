#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "regression"))
from gcta_runner import run_gcta
from plink import write_plink
from synthetic import family_markers, write_family_phenotypes


def write_synthetic_plink(prefix, sample_count=2048, marker_count=8193, profile="regular"):
	prefix = Path(prefix)
	prefix.parent.mkdir(parents=True, exist_ok=True)
	prefix.with_suffix(".fam").write_text(
		"".join(f"F1 I{sample} 0 0 0 -9\n" for sample in range(sample_count)), encoding="utf-8"
	)
	prefix.with_suffix(".bim").write_text(
		"".join(f"1 rs{marker} 0 {1000 + marker} A G\n" for marker in range(marker_count)),
		encoding="utf-8",
	)
	codes = {0: 3, 1: 2, 2: 0, None: 1}
	missing_rates = (0, 10, 50, 200)
	with prefix.with_suffix(".bed").open("wb") as stream:
		stream.write(b"\x6c\x1b\x01")
		for marker in range(marker_count):
			for offset in range(0, sample_count, 4):
				packed = 0
				for within_byte, sample in enumerate(range(offset, min(offset + 4, sample_count))):
					genotype = (sample * 3 + marker * 5 + sample * marker) % 3
					if profile == "missing":
						if marker % 257 == 0:
							genotype = marker % 3
						elif marker % 131 == 0:
							genotype = 1 if sample == 0 else 0
						missing_value = (sample * 31 + marker * 17 + sample * marker) % 1000
						if missing_value < missing_rates[marker % 4]:
							genotype = None
					packed |= codes[genotype] << (2 * within_byte)
				stream.write(bytes((packed,)))


def prepare_greml(gcta, data, work):
	legacy = data / "gcta_legacy" / "test"
	keep = work / "input" / "keep1024.txt"
	keep.parent.mkdir(parents=True, exist_ok=True)
	with legacy.with_suffix(".fam").open(encoding="utf-8") as stream:
		keep.write_text("".join(stream.readlines()[:1024]), encoding="utf-8")
	grm = run_gcta(
		gcta,
		work / "setup",
		"grm",
		["--bfile", legacy, "--keep", keep, "--make-grm"],
	)
	return [
		"--grm", grm,
		"--pheno", legacy.with_suffix(".phen"),
		"--reml",
		"--reml-priors", 0.5,
		"--reml-maxit", 100,
		"--reml-no-lrt",
		"--reml-est-fix",
	], 1024, 1000


def prepare_fastgwa(gcta, work, scenario):
	sample_count = 2048 if scenario.endswith("_large") else 1024
	marker_count = 16385 if scenario.endswith("_large") else 8193
	case_fraction = 0.1 if scenario == "fastgwa_glmm_imbalanced" else 0.5
	markers, chromosomes = family_markers(sample_count=sample_count, marker_count=marker_count)
	prefix = work / "input" / "family"
	write_plink(prefix, markers, chromosomes)
	pheno, qcovar = write_family_phenotypes(prefix, markers, case_fraction=case_fraction)
	arguments = ["--bfile", prefix]
	if scenario not in ("fastgwa_lr", "fastgwa_lr_large"):
		sparse = run_gcta(
			gcta,
			work / "setup",
			"sparse",
			["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05],
		)
		arguments.extend(("--grm-sparse", sparse))
	option = {
		"fastgwa_lr": "--fastGWA-lr",
		"fastgwa_mlm": "--fastGWA-mlm",
		"fastgwa_exact": "--fastGWA-mlm-exact",
		"fastgwa_glmm": "--fastGWA-mlm-binary",
		"fastgwa_lr_large": "--fastGWA-lr",
		"fastgwa_mlm_large": "--fastGWA-mlm",
		"fastgwa_glmm_imbalanced": "--fastGWA-mlm-binary",
	}[scenario]
	arguments.extend((option, "--pheno", pheno, "--qcovar", qcovar, "--seed", 12345))
	if scenario in ("fastgwa_glmm", "fastgwa_glmm_imbalanced"):
		arguments.extend(("--mpheno", 2))
	return arguments, sample_count, marker_count


def prepare_scenario(gcta, data, work, scenario):
	if scenario == "grm_legacy":
		prefix = data / "gcta_legacy" / "test"
		return ["--bfile", prefix, "--make-grm"], 3925, 1000
	if scenario == "grm_synthetic":
		prefix = work / "input" / "synthetic"
		write_synthetic_plink(prefix)
		return ["--bfile", prefix, "--make-grm"], 2048, 8193
	if scenario == "grm_large":
		prefix = work / "input" / "synthetic"
		write_synthetic_plink(prefix, sample_count=4096, marker_count=16385)
		return ["--bfile", prefix, "--make-grm"], 4096, 16385
	if scenario == "grm_missing":
		prefix = work / "input" / "synthetic"
		write_synthetic_plink(prefix, profile="missing")
		return ["--bfile", prefix, "--make-grm"], 2048, 8193
	if scenario == "greml_legacy":
		return prepare_greml(gcta, data, work)
	if scenario in ("fastgwa_lr_binary", "fastgwa_glmm_binary"):
		base_scenario = "fastgwa_lr" if scenario == "fastgwa_lr_binary" else "fastgwa_glmm"
		arguments, sample_count, marker_count = prepare_fastgwa(gcta, work, base_scenario)
		arguments.append("--save-bin")
		return arguments, sample_count, marker_count
	return prepare_fastgwa(gcta, work, scenario)


def sha256(path):
	digest = hashlib.sha256()
	with Path(path).open("rb") as stream:
		for block in iter(lambda: stream.read(1024 * 1024), b""):
			digest.update(block)
	return digest.hexdigest()


def cpu_model():
	for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
		if line.startswith("model name"):
			return line.split(":", 1)[1].strip()
	return "unknown"


def source_state():
	root = Path(__file__).resolve().parents[2]
	commit = subprocess.run(
		["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True, check=True
	).stdout.strip()
	dirty = subprocess.run(
		["git", "-C", root, "status", "--porcelain"], capture_output=True, text=True, check=True
	).stdout.strip()
	return commit, bool(dirty)


def measure(command, work_dir, name, threads):
	work_dir.mkdir(parents=True, exist_ok=True)
	console_path = work_dir / f"{name}.console"
	environment = os.environ.copy()
	environment["OMP_NUM_THREADS"] = str(threads)
	environment["OMP_DYNAMIC"] = "FALSE"
	environment["MKL_NUM_THREADS"] = str(threads)
	environment["MKL_DYNAMIC"] = "FALSE"
	with console_path.open("w", encoding="utf-8") as console:
		start = time.perf_counter()
		process = subprocess.Popen(
			[str(argument) for argument in command],
			cwd=work_dir,
			env=environment,
			stdout=console,
			stderr=subprocess.STDOUT,
		)
		_, status, usage = os.wait4(process.pid, 0)
		elapsed = time.perf_counter() - start
	process.returncode = os.waitstatus_to_exitcode(status)
	assert process.returncode == 0, f"Benchmark failed; see {console_path}"
	return {
		"run": name,
		"elapsed_seconds": elapsed,
		"user_cpu_seconds": usage.ru_utime,
		"system_cpu_seconds": usage.ru_stime,
		"max_rss_kb": usage.ru_maxrss,
	}


def main():
	parser = argparse.ArgumentParser(description="Run reproducible GCTA performance benchmarks")
	parser.add_argument("--gcta", required=True)
	parser.add_argument("--data", required=True)
	parser.add_argument("--work", required=True)
	parser.add_argument(
		"--scenario",
		required=True,
		choices=(
			"grm_legacy",
			"grm_synthetic",
			"grm_large",
			"grm_missing",
			"greml_legacy",
			"fastgwa_lr",
			"fastgwa_mlm",
			"fastgwa_exact",
			"fastgwa_glmm",
			"fastgwa_lr_large",
			"fastgwa_mlm_large",
			"fastgwa_glmm_imbalanced",
			"fastgwa_lr_binary",
			"fastgwa_glmm_binary",
		),
	)
	parser.add_argument("--threads", nargs="+", type=int, default=(1, 2, 4))
	parser.add_argument("--repeat", type=int, default=5)
	args = parser.parse_args()
	assert args.repeat > 0 and all(thread > 0 for thread in args.threads)
	gcta = Path(args.gcta).resolve()
	work = Path(args.work).resolve()
	data = Path(args.data).resolve()
	scenario_work = work / args.scenario
	arguments, sample_count, marker_count = prepare_scenario(
		gcta, data, scenario_work, args.scenario
	)

	source_commit, source_dirty = source_state()
	results = {
		"scenario": args.scenario,
		"source_commit": source_commit,
		"source_dirty": source_dirty,
		"sample_count": sample_count,
		"marker_count": marker_count,
		"repeat": args.repeat,
		"binary": str(gcta),
		"binary_sha256": sha256(gcta),
		"platform": platform.platform(),
		"cpu": cpu_model(),
		"python": platform.python_version(),
		"measurements": {},
	}
	for threads in args.threads:
		thread_work = scenario_work / f"threads_{threads}"
		output = thread_work / "benchmark"
		command = [
			gcta,
			*arguments,
			"--out", output,
			"--thread-num", threads,
		]
		measure(command, thread_work, "warmup", threads)
		measurements = [
			measure(command, thread_work, f"run_{run}", threads) for run in range(args.repeat)
		]
		results["measurements"][str(threads)] = measurements
		elapsed = [measurement["elapsed_seconds"] for measurement in measurements]
		print(
			f"{args.scenario}, threads={threads}: median={statistics.median(elapsed):.6f}s "
			f"range={min(elapsed):.6f}-{max(elapsed):.6f}s, "
			f"peak RSS={max(item['max_rss_kb'] for item in measurements)} KiB"
		)
	results_path = work / f"{args.scenario}.json"
	results_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
	print(f"Saved raw benchmark measurements to {results_path}")


if __name__ == "__main__":
	main()
