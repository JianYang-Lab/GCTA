#!/usr/bin/env python3

import argparse
import json
import platform
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "regression"))
from run import cpu_model, measure, prepare_scenario, sha256, source_state

DEFAULT_SCENARIOS = (
    "grm_legacy",
    "grm_synthetic",
    "grm_missing",
    "greml_legacy",
    "fastgwa_lr",
    "fastgwa_mlm",
    "fastgwa_exact",
    "fastgwa_glmm",
    "fastgwa_2df",
    "fastgwa_2df_nosandwich",
    "fastgwa_2df_mixed",
    "fastgwa_2df_mixed_nosandwich",
    "fastgwa_glmm_large",
    "ld_matrix",
)


def parse_args():
    parser = argparse.ArgumentParser(description="Compare GCTA benchmark binaries")
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--data", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenarios", nargs="+", default=DEFAULT_SCENARIOS)
    parser.add_argument("--threads", nargs="+", type=int, default=(1, 4))
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--max-slowdown", type=float, default=0.10)
    parser.add_argument("--max-rss-increase", type=float, default=0.10)
    parser.add_argument("--min-wall-seconds", type=float, default=0.25)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--cpu-affinity", default=None)
    args = parser.parse_args()
    assert args.repeat > 0 and all(thread > 0 for thread in args.threads)
    if args.cpu_affinity is not None:
        args.cpu_affinity = [int(cpu) for cpu in args.cpu_affinity.split(",") if cpu]
        assert len(args.cpu_affinity) >= max(args.threads)
    return args


def run_once(binary, arguments, work_dir, name, threads, cpus=None):
    output = work_dir / "benchmark"
    command = [
        binary,
        *arguments,
        "--out", output,
        "--thread-num", threads,
    ]
    return measure(command, work_dir, name, threads, cpus)


def compare_scenario(args, baseline, candidate, data, work, scenario):
    scenario_work = work / scenario
    arguments, sample_count, marker_count = prepare_scenario(
        baseline, data, scenario_work, scenario
    )
    results = []
    for threads in args.threads:
        baseline_work = scenario_work / "baseline" / f"threads_{threads}"
        candidate_work = scenario_work / "candidate" / f"threads_{threads}"
        baseline_work.mkdir(parents=True, exist_ok=True)
        candidate_work.mkdir(parents=True, exist_ok=True)
        run_once(baseline, arguments, baseline_work, "warmup", threads, args.cpu_affinity)
        run_once(candidate, arguments, candidate_work, "warmup", threads, args.cpu_affinity)

        baseline_runs = []
        candidate_runs = []
        for run in range(args.repeat):
            if run % 2 == 0:
                baseline_runs.append(
                    run_once(baseline, arguments, baseline_work, f"run_{run}", threads, args.cpu_affinity)
                )
                candidate_runs.append(
                    run_once(candidate, arguments, candidate_work, f"run_{run}", threads, args.cpu_affinity)
                )
            else:
                candidate_runs.append(
                    run_once(candidate, arguments, candidate_work, f"run_{run}", threads, args.cpu_affinity)
                )
                baseline_runs.append(
                    run_once(baseline, arguments, baseline_work, f"run_{run}", threads, args.cpu_affinity)
                )

        baseline_time = statistics.median(item["elapsed_seconds"] for item in baseline_runs)
        candidate_time = statistics.median(item["elapsed_seconds"] for item in candidate_runs)
        baseline_cpu = statistics.median(
            item["user_cpu_seconds"] + item["system_cpu_seconds"] for item in baseline_runs
        )
        candidate_cpu = statistics.median(
            item["user_cpu_seconds"] + item["system_cpu_seconds"] for item in candidate_runs
        )
        baseline_rss = statistics.median(item["max_rss_kb"] for item in baseline_runs)
        candidate_rss = statistics.median(item["max_rss_kb"] for item in candidate_runs)
        time_ratio = candidate_time / baseline_time if baseline_time > 0 else 1.0
        cpu_ratio = candidate_cpu / baseline_cpu if baseline_cpu > 0 else 1.0
        rss_ratio = candidate_rss / baseline_rss if baseline_rss > 0 else 1.0
        wall_regression = (
            baseline_time >= args.min_wall_seconds
            and time_ratio > 1.0 + args.max_slowdown
        )
        failed = (
            cpu_ratio > 1.0 + args.max_slowdown
            or wall_regression
            or rss_ratio > 1.0 + args.max_rss_increase
        )
        results.append({
            "threads": threads,
            "baseline_median_seconds": baseline_time,
            "candidate_median_seconds": candidate_time,
            "time_ratio": time_ratio,
            "baseline_median_cpu_seconds": baseline_cpu,
            "candidate_median_cpu_seconds": candidate_cpu,
            "cpu_ratio": cpu_ratio,
            "baseline_median_rss_kb": baseline_rss,
            "candidate_median_rss_kb": candidate_rss,
            "rss_ratio": rss_ratio,
            "status": "fail" if failed else "pass",
        })
        print(
            f"{scenario}, threads={threads}: "
            f"wall {baseline_time:.6f}s -> {candidate_time:.6f}s "
            f"(x{time_ratio:.3f}), CPU {baseline_cpu:.6f}s -> {candidate_cpu:.6f}s "
            f"(x{cpu_ratio:.3f}), RSS {baseline_rss:.0f} -> {candidate_rss:.0f} KiB "
            f"(x{rss_ratio:.3f}) [{results[-1]['status']}]"
        )
    return {
        "scenario": scenario,
        "sample_count": sample_count,
        "marker_count": marker_count,
        "measurements": results,
    }


def main():
    args = parse_args()
    baseline = Path(args.baseline).resolve()
    candidate = Path(args.candidate).resolve()
    data = Path(args.data).resolve()
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    commit, dirty = source_state()
    report = {
        "baseline": str(baseline),
        "baseline_sha256": sha256(baseline),
        "candidate": str(candidate),
        "candidate_sha256": sha256(candidate),
        "source_commit": commit,
        "source_dirty": dirty,
        "platform": platform.platform(),
        "cpu": cpu_model(),
        "repeat": args.repeat,
        "max_slowdown": args.max_slowdown,
        "max_rss_increase": args.max_rss_increase,
        "min_wall_seconds": args.min_wall_seconds,
        "cpu_affinity": args.cpu_affinity,
        "results": [],
    }
    for scenario in args.scenarios:
        report["results"].append(
            compare_scenario(args, baseline, candidate, data, work, scenario)
        )
    report_path = work / "compare.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    failed = any(
        measurement["status"] == "fail"
        for result in report["results"]
        for measurement in result["measurements"]
    )
    print(f"Saved benchmark comparison to {report_path}")
    if failed and not args.report_only:
        print("Benchmark regression detected.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
