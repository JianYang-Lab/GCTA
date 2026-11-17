#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, assert_tables_close, read_floats, run_gcta
from plink import expected_grm, synthetic_markers, write_plink
from synthetic import family_markers


def make_fixture(work_dir):
	prefix = Path(work_dir) / "input" / "synthetic"
	sample_ids = write_plink(prefix, synthetic_markers())
	return prefix, sample_ids


def check_oracle(gcta, work_dir, algorithm, markers=None):
	if markers is None:
		markers = synthetic_markers()
	prefix = Path(work_dir) / "input" / "synthetic"
	sample_ids = write_plink(prefix, markers)
	output = run_gcta(
		gcta,
		work_dir,
		f"algorithm_{algorithm}",
		["--bfile", prefix, "--make-grm", "--make-grm-alg", algorithm],
	)
	expected_grm_values, expected_counts = expected_grm(markers, algorithm)
	assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected_grm_values, output)
	actual_counts = read_floats(output.with_suffix(".grm.N.bin"))
	assert actual_counts == expected_counts, f"GRM marker counts differ: {output}"
	expected_ids = "".join(f"{fid}\t{iid}\n" for fid, iid in sample_ids)
	assert output.with_suffix(".grm.id").read_text(encoding="utf-8") == expected_ids


def check_missing(gcta, work_dir):
	markers = [
		[0, 1, None],
		[2, 1, 0],
		[0, None, 2],
		[None, 2, 1],
		[0, 0, 0],
		[None, None, 1],
	]
	for algorithm in (0, 1):
		check_oracle(gcta, Path(work_dir) / f"algorithm_{algorithm}", algorithm, markers)


def check_threads(gcta, data_dir, work_dir):
	prefix = Path(data_dir) / "gcta_legacy" / "test"
	single = run_gcta(gcta, work_dir, "single", ["--bfile", prefix, "--make-grm"], threads=1)
	multi = run_gcta(gcta, work_dir, "multi", ["--bfile", prefix, "--make-grm"], threads=4)
	assert_grm_equal(single, multi)


def check_sparse(gcta, work_dir):
	markers, chromosomes = family_markers()
	prefix = Path(work_dir) / "input" / "family"
	write_plink(prefix, markers, chromosomes)
	dense = run_gcta(gcta, work_dir, "dense", ["--bfile", prefix, "--make-grm"])
	direct = run_gcta(
		gcta, work_dir, "direct", ["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05]
	)
	multi = run_gcta(
		gcta,
		work_dir,
		"multi",
		["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05],
		threads=4,
	)
	pruned = run_gcta(
		gcta, work_dir, "pruned", ["--grm", dense, "--make-bK-sparse", 0.05]
	)
	for suffix in ("grm.sp", "grm.id"):
		assert_tables_close(direct.with_suffix(f".{suffix}"), pruned.with_suffix(f".{suffix}"))
		assert_tables_close(direct.with_suffix(f".{suffix}"), multi.with_suffix(f".{suffix}"))
	rows = [line.split() for line in direct.with_suffix(".grm.sp").read_text().splitlines()]
	assert len(rows) > len(markers[0])
	for row, column, value in rows:
		row = int(row)
		column = int(column)
		assert 0 <= column <= row < len(markers[0])
		assert float(value) >= 0.05


def check_snp_merge(gcta, work_dir):
	prefix, _ = make_fixture(work_dir)
	whole = run_gcta(gcta, work_dir, "whole", ["--bfile", prefix, "--make-grm"])
	subset_outputs = []
	for subset, remainder in enumerate((0, 1)):
		snplist = Path(work_dir) / f"subset_{subset}.snplist"
		snplist.write_text(
			"".join(f"rs{index}\n" for index in range(131) if index % 2 == remainder), encoding="utf-8"
		)
		subset_outputs.append(run_gcta(
			gcta,
			work_dir,
			f"subset_{subset}",
			["--bfile", prefix, "--extract", snplist, "--make-grm"],
		))
	mgrm = Path(work_dir) / "mgrm.txt"
	mgrm.write_text("".join(f"{output}\n" for output in subset_outputs), encoding="utf-8")
	merged = run_gcta(gcta, work_dir, "merged", ["--mgrm", mgrm, "--make-grm"])
	assert_grm_equal(whole, merged)


def check_partitions(gcta, work_dir):
	prefix, _ = make_fixture(work_dir)
	whole = run_gcta(gcta, work_dir, "whole", ["--bfile", prefix, "--make-grm"])
	part_outputs = []
	for part in range(1, 4):
		base = run_gcta(
			gcta,
			work_dir,
			f"partition_{part}",
			["--bfile", prefix, "--make-grm-part", 3, part],
			threads=2,
		)
		part_outputs.append(Path(f"{base}.part_3_{part}"))
	joined = Path(work_dir) / "joined"
	for suffix in ("grm.bin", "grm.N.bin", "grm.id"):
		Path(f"{joined}.{suffix}").write_bytes(
			b"".join(Path(f"{part}.{suffix}").read_bytes() for part in part_outputs)
		)
	assert_grm_equal(whole, joined)


def main():
	parser = argparse.ArgumentParser(description="GCTA GRM regression tests")
	parser.add_argument("--gcta", required=True)
	parser.add_argument("--data", required=True)
	parser.add_argument("--work", required=True)
	parser.add_argument(
		"--scenario",
		required=True,
		choices=("standard", "weighted", "missing", "threads", "sparse", "snp_merge", "partitions"),
	)
	args = parser.parse_args()
	work_dir = Path(args.work)
	work_dir.mkdir(parents=True, exist_ok=True)

	if args.scenario == "standard":
		check_oracle(args.gcta, work_dir, 0)
	elif args.scenario == "weighted":
		check_oracle(args.gcta, work_dir, 1)
	elif args.scenario == "missing":
		check_missing(args.gcta, work_dir)
	elif args.scenario == "threads":
		check_threads(args.gcta, args.data, work_dir)
	elif args.scenario == "sparse":
		check_sparse(args.gcta, work_dir)
	elif args.scenario == "snp_merge":
		check_snp_merge(args.gcta, work_dir)
	else:
		check_partitions(args.gcta, work_dir)
	print(f"Passed GRM regression scenario: {args.scenario}")


if __name__ == "__main__":
	main()
