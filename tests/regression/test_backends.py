#!/usr/bin/env python3

import argparse
import ctypes.util
import os
from pathlib import Path
import platform

from gcta_runner import assert_floats_close, read_floats, run_gcta
from plink import expected_grm, synthetic_markers, write_plink
from test_errors import expect_error


def main():
	parser = argparse.ArgumentParser(description="GCTA numerical backend regression test")
	parser.add_argument("--gcta", required=True)
	parser.add_argument("--data", required=True)
	parser.add_argument("--work", required=True)
	parser.add_argument("--scenario", default="grm", choices=("grm",))
	args = parser.parse_args()
	if platform.system() != "Linux":
		print("SKIP: runtime backend selection is currently implemented on Linux")
		return
	markers = synthetic_markers()
	for backend in ("mkl", "openblas"):
		if backend == "openblas" and not ctypes.util.find_library("openblaso"):
			print("SKIP: OpenBLAS is unavailable")
			continue
		for threads in (1, 4):
			for algorithm in (0, 1):
				work = Path(args.work) / backend / f"threads_{threads}" / f"algorithm_{algorithm}"
				prefix = work / "input" / "synthetic"
				sample_ids = write_plink(prefix, markers)
				os.environ["GCTA_NUMERIC_BACKEND"] = backend
				output = run_gcta(
					args.gcta,
					work,
					"grm",
					["--bfile", prefix, "--make-grm", "--make-grm-alg", algorithm],
					threads=threads,
				)
				expected, counts = expected_grm(markers, algorithm)
				assert_floats_close(read_floats(Path(f"{output}.grm.bin")), expected, output)
				assert read_floats(Path(f"{output}.grm.N.bin")) == counts
				expected_ids = "".join(f"{fid}\t{iid}\n" for fid, iid in sample_ids)
				assert Path(f"{output}.grm.id").read_text(encoding="utf-8") == expected_ids
				log = Path(f"{output}.log").read_text(encoding="utf-8")
				assert f"GRM BLAS backend: {backend.upper()}" in log
				assert "Error:" not in log
	os.environ["GCTA_NUMERIC_BACKEND"] = "invalid"
	invalid_work = Path(args.work) / "invalid"
	expect_error(
		args.gcta,
		invalid_work,
		"invalid_backend",
		["--bfile", prefix, "--make-grm", "--out", invalid_work / "result"],
	)
	print("Passed numerical backend regression scenario: grm")


if __name__ == "__main__":
	main()
