#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "regression"))
from gcta_runner import run_gcta
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


def write_plink_binary(prefix, genotypes, chromosomes=None):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    marker_count = len(genotypes)
    sample_count = len(genotypes[0])
    if chromosomes is None:
        chromosomes = [1] * marker_count
    prefix.with_suffix(".fam").write_text(
        "".join(f"F1 I{index} 0 0 0 -9\n" for index in range(sample_count)),
        encoding="utf-8",
    )
    prefix.with_suffix(".bim").write_text(
        "".join(
            f"{chromosomes[index]} rs{index} 0 {1000 + index} A G\n"
            for index in range(marker_count)
        ),
        encoding="utf-8",
    )
    code_by_genotype = {2: 0, None: 1, 1: 2, 0: 3}
    with prefix.with_suffix(".bed").open("wb") as stream:
        stream.write(b"\x6c\x1b\x01")
        for marker in genotypes:
            for offset in range(0, sample_count, 4):
                value = 0
                for within_byte, genotype in enumerate(marker[offset:offset + 4]):
                    value |= code_by_genotype[genotype] << (2 * within_byte)
                stream.write(bytes((value,)))


def repeat_plink(source, destination, repeats):
    source = Path(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    bim_lines = source.with_suffix(".bim").read_text(encoding="utf-8").splitlines()
    marker_count = len(bim_lines)
    bed = source.with_suffix(".bed").read_bytes()
    assert (len(bed) - 3) % marker_count == 0, f"Unexpected BED size: {source}"
    marker_bytes = (len(bed) - 3) // marker_count
    destination.with_suffix(".bed").write_bytes(bed[:3] + bed[3:] * repeats)
    destination.with_suffix(".fam").write_bytes(source.with_suffix(".fam").read_bytes())
    lines = []
    for repeat in range(repeats):
        for index, line in enumerate(bim_lines):
            fields = line.split()
            fields[1] = f"rs{repeat * marker_count + index}"
            fields[3] = str(1000 + index + repeat * 10000000)
            lines.append(" ".join(fields))
    destination.with_suffix(".bim").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_benchmark_phenotypes(prefix, sample_count, seed=20261119, binary=True):
    prefix = Path(prefix)
    rng = random.Random(seed)
    quantitative = [rng.gauss(0.0, 1.0) for _ in range(sample_count)]
    binary_values = [1 if rng.random() < 0.5 else 0 for _ in range(sample_count)]
    with prefix.with_suffix(".phen").open("w", encoding="utf-8") as stream:
        for sample in range(sample_count):
            stream.write(
                f"F1 I{sample}\t{quantitative[sample]:.17g}\t{binary_values[sample]}\n"
            )
    with prefix.with_suffix(".qcovar").open("w", encoding="utf-8") as stream:
        for sample in range(sample_count):
            stream.write(
                f"F1 I{sample}\t{rng.gauss(0.0, 1.0):.17g}\t{rng.gauss(0.0, 1.0):.17g}\n"
            )
    envir = prefix.parent / "envir.txt"
    with envir.open("w", encoding="utf-8") as stream:
        for sample in range(sample_count):
            stream.write(f"F1 I{sample}\t{rng.gauss(0.0, 1.0):.17g}\n")
    return prefix.with_suffix(".phen"), prefix.with_suffix(".qcovar"), envir


def prepare_fastgwa_2df(gcta, work, sandwich=True, sample_count=4096, marker_count=8193, repeats=4):
    base = work / "input" / "base"
    write_synthetic_plink(base, sample_count=sample_count, marker_count=marker_count)
    prefix = work / "input" / "synthetic"
    repeat_plink(base, prefix, repeats)
    pheno, qcovar, envir = write_benchmark_phenotypes(prefix, sample_count)
    arguments = [
        "--bfile", prefix,
        "--fastGWA-lr",
        "--pheno", pheno,
        "--qcovar", qcovar,
        "--envir", envir,
        "--seed", 12345,
    ]
    if not sandwich:
        arguments.append("--noSandwich")
    return arguments, sample_count, marker_count * repeats


def prepare_fastgwa_family(gcta, work, option, envir, sample_count=1024, marker_count=8193, repeats=4, extra=()):
    markers, chromosomes = family_markers(sample_count=sample_count, marker_count=marker_count)
    base = work / "input" / "base"
    write_plink_binary(base, markers, chromosomes)
    prefix = work / "input" / "family"
    repeat_plink(base, prefix, repeats)
    pheno, qcovar, envir_path = write_benchmark_phenotypes(prefix, sample_count)
    sparse = run_gcta(
        gcta,
        work / "setup",
        "sparse",
        ["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05],
    )
    arguments = [
        "--bfile", prefix,
        "--grm-sparse", sparse,
        option,
        "--pheno", pheno,
        "--qcovar", qcovar,
        "--seed", 12345,
    ]
    if envir:
        arguments.extend(("--envir", envir_path))
    arguments.extend(extra)
    return arguments, sample_count, marker_count * repeats


def prepare_ld_matrix(work, sample_count=1024, marker_count=8193):
    prefix = work / "input" / "synthetic"
    write_synthetic_plink(prefix, sample_count=sample_count, marker_count=marker_count)
    lines = prefix.with_suffix(".bim").read_text(encoding="utf-8").splitlines()
    spaced = []
    for index, line in enumerate(lines):
        fields = line.split()
        fields[3] = str(1000 + index * 1000)
        spaced.append(" ".join(fields))
    prefix.with_suffix(".bim").write_text("\n".join(spaced) + "\n", encoding="utf-8")
    return ["--bfile", prefix, "--ld-matrix", "--ld-wind", "1000"], sample_count, marker_count


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
    write_plink_binary(prefix, markers, chromosomes)
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
    if scenario in ("fastgwa_2df", "fastgwa_2df_nosandwich"):
        return prepare_fastgwa_2df(gcta, work, sandwich=scenario == "fastgwa_2df")
    if scenario in ("fastgwa_2df_mixed", "fastgwa_2df_mixed_nosandwich"):
        extra = () if scenario == "fastgwa_2df_mixed" else ("--noSandwich",)
        return prepare_fastgwa_family(gcta, work, "--fastGWA-mlm", True, extra=extra)
    if scenario == "fastgwa_glmm_large":
        return prepare_fastgwa_family(
            gcta, work, "--fastGWA-mlm-binary", False, extra=("--mpheno", "2")
        )
    if scenario == "ld_matrix":
        return prepare_ld_matrix(work)
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


# Running the benchmark in this process would make ru_maxrss inherit the
# Python harness high-water mark (fork/exec keeps it), so a small helper does
# the fork and reports the child's own resource usage.
MEASURE_HELPER = r"""
import json
import os
import sys
import time

result_path = sys.argv[1]
command = json.loads(sys.argv[2])
start = time.perf_counter()
pid = os.fork()
if pid == 0:
    try:
        os.execvpe(command[0], command, os.environ)
    except Exception:
        os._exit(127)
_, status, usage = os.wait4(pid, 0)
elapsed = time.perf_counter() - start
with open(result_path, "w", encoding="utf-8") as stream:
    json.dump({
        "returncode": os.waitstatus_to_exitcode(status),
        "elapsed_seconds": elapsed,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "max_rss_kb": usage.ru_maxrss,
    }, stream)
"""


def measure(command, work_dir, name, threads, cpus=None):
    work_dir.mkdir(parents=True, exist_ok=True)
    console_path = work_dir / f"{name}.console"
    result_path = work_dir / f"{name}.measure.json"
    if cpus:
        command = ["taskset", "-c", ",".join(str(cpu) for cpu in cpus[:threads]), *command]
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(threads)
    environment["OMP_DYNAMIC"] = "FALSE"
    environment["MKL_NUM_THREADS"] = str(threads)
    environment["MKL_DYNAMIC"] = "FALSE"
    with console_path.open("w", encoding="utf-8") as console:
        subprocess.run(
            [
                sys.executable, "-c", MEASURE_HELPER,
                str(result_path), json.dumps([str(argument) for argument in command]),
            ],
            cwd=work_dir,
            env=environment,
            stdout=console,
            stderr=subprocess.STDOUT,
            check=True,
        )
    measured = json.loads(result_path.read_text(encoding="utf-8"))
    assert measured["returncode"] == 0, f"Benchmark failed; see {console_path}"
    # Keep only the console log; some scenarios (e.g. LD matrices) produce
    # hundreds of megabytes per run and the comparison only needs timings.
    for entry in work_dir.iterdir():
        if entry != console_path and entry.is_file():
            entry.unlink()
    return {
        "run": name,
        "elapsed_seconds": measured["elapsed_seconds"],
        "user_cpu_seconds": measured["user_cpu_seconds"],
        "system_cpu_seconds": measured["system_cpu_seconds"],
        "max_rss_kb": measured["max_rss_kb"],
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
            "fastgwa_2df",
            "fastgwa_2df_nosandwich",
            "fastgwa_2df_mixed",
            "fastgwa_2df_mixed_nosandwich",
            "fastgwa_glmm_large",
            "ld_matrix",
        ),
    )
    parser.add_argument("--threads", nargs="+", type=int, default=(1, 2, 4))
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--cpu-affinity", default=None)
    args = parser.parse_args()
    assert args.repeat > 0 and all(thread > 0 for thread in args.threads)
    cpus = None
    if args.cpu_affinity is not None:
        cpus = [int(cpu) for cpu in args.cpu_affinity.split(",") if cpu]
        assert len(cpus) >= max(args.threads)
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
        "cpu_affinity": args.cpu_affinity,
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
        measure(command, thread_work, "warmup", threads, cpus)
        measurements = [
            measure(command, thread_work, f"run_{run}", threads, cpus) for run in range(args.repeat)
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
