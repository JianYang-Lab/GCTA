#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time


def write_synthetic_plink(prefix, sample_count=2048, marker_count=8193):
	prefix = Path(prefix)
	prefix.parent.mkdir(parents=True, exist_ok=True)
	prefix.with_suffix(".fam").write_text(
		"".join(f"F1 I{sample} 0 0 0 -9\n" for sample in range(sample_count)), encoding="utf-8"
	)
	prefix.with_suffix(".bim").write_text(
		"".join(f"1 rs{marker} 0 {1000 + marker} A G\n" for marker in range(marker_count)),
		encoding="utf-8",
	)
	codes = (3, 2, 0)
	with prefix.with_suffix(".bed").open("wb") as stream:
		stream.write(b"\x6c\x1b\x01")
		for marker in range(marker_count):
			for offset in range(0, sample_count, 4):
				packed = 0
				for within_byte, sample in enumerate(range(offset, min(offset + 4, sample_count))):
					genotype = (sample * 3 + marker * 5 + sample * marker) % 3
					packed |= codes[genotype] << (2 * within_byte)
				stream.write(bytes((packed,)))


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
	parser.add_argument("--scenario", required=True, choices=("grm_legacy", "grm_synthetic"))
	parser.add_argument("--threads", nargs="+", type=int, default=(1, 2, 4))
	parser.add_argument("--repeat", type=int, default=5)
	args = parser.parse_args()
	assert args.repeat > 0 and all(thread > 0 for thread in args.threads)
	gcta = Path(args.gcta).resolve()
	work = Path(args.work).resolve()
	data = Path(args.data).resolve()
	if args.scenario == "grm_legacy":
		prefix = data / "gcta_legacy" / "test"
		sample_count, marker_count = 3925, 1000
	else:
		prefix = work / "input" / "synthetic"
		write_synthetic_plink(prefix)
		sample_count, marker_count = 2048, 8193

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
		thread_work = work / f"threads_{threads}"
		output = thread_work / "benchmark"
		command = [
			gcta,
			"--bfile", prefix,
			"--make-grm",
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
