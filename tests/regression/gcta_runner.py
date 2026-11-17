import math
import os
from pathlib import Path
import struct
import subprocess


def run_gcta(gcta, work_dir, name, arguments, threads=1, timeout=120):
	work_dir = Path(work_dir).resolve()
	work_dir.mkdir(parents=True, exist_ok=True)
	output = work_dir / name
	environment = os.environ.copy()
	environment["OMP_NUM_THREADS"] = str(threads)
	environment["OMP_DYNAMIC"] = "FALSE"
	environment["MKL_NUM_THREADS"] = str(threads)
	environment["MKL_DYNAMIC"] = "FALSE"
	command = [
		str(Path(gcta).resolve()),
		*(str(argument.resolve() if isinstance(argument, Path) else argument) for argument in arguments),
		"--out", str(output),
		"--thread-num", str(threads),
	]
	result = subprocess.run(
		command,
		cwd=work_dir,
		env=environment,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT,
		text=True,
		timeout=timeout,
		check=False,
	)
	(work_dir / f"{name}.console").write_text(result.stdout, encoding="utf-8")
	assert result.returncode == 0, (
		f"GCTA failed with exit code {result.returncode}: {' '.join(command)}\n{result.stdout}"
	)
	return output


def read_floats(path):
	data = Path(path).read_bytes()
	assert len(data) % 4 == 0, f"Invalid float file size: {path}"
	return [item[0] for item in struct.iter_unpack("<f", data)]


def assert_floats_close(actual, expected, path, rel_tol=2e-6, abs_tol=2e-6):
	assert len(actual) == len(expected), (
		f"Element count mismatch for {path}: {len(actual)} != {len(expected)}"
	)
	for index, (observed, reference) in enumerate(zip(actual, expected)):
		assert math.isfinite(observed), f"Non-finite value at {path}[{index}]: {observed}"
		assert math.isclose(observed, reference, rel_tol=rel_tol, abs_tol=abs_tol), (
			f"Value mismatch at {path}[{index}]: {observed} != {reference}"
		)


def read_table(path):
	return [line.split() for line in Path(path).read_text(encoding="utf-8").splitlines()]


def assert_tables_close(actual, expected, rel_tol=2e-5, abs_tol=1e-6):
	actual_rows = read_table(actual)
	expected_rows = read_table(expected)
	assert len(actual_rows) == len(expected_rows), (
		f"Row count mismatch for {actual}: {len(actual_rows)} != {len(expected_rows)}"
	)
	for row_index, (actual_row, expected_row) in enumerate(zip(actual_rows, expected_rows)):
		assert len(actual_row) == len(expected_row), f"Column count mismatch at {actual}:{row_index + 1}"
		for column_index, (observed, reference) in enumerate(zip(actual_row, expected_row)):
			if observed == reference:
				continue
			try:
				observed_number = float(observed)
				reference_number = float(reference)
			except ValueError:
				raise AssertionError(
					f"Text mismatch at {actual}:{row_index + 1}:{column_index + 1}: "
					f"{observed} != {reference}"
				)
			if math.isnan(observed_number) and math.isnan(reference_number):
				continue
			assert math.isclose(observed_number, reference_number, rel_tol=rel_tol, abs_tol=abs_tol), (
				f"Value mismatch at {actual}:{row_index + 1}:{column_index + 1}: "
				f"{observed} != {reference}"
			)


def assert_grm_equal(first, second, rel_tol=2e-6, abs_tol=2e-6):
	first = Path(first)
	second = Path(second)
	for suffix in ("grm.bin", "grm.N.bin"):
		first_path = Path(f"{first}.{suffix}")
		second_path = Path(f"{second}.{suffix}")
		first_data = first_path.read_bytes()
		second_data = second_path.read_bytes()
		if first_data == second_data:
			continue
		assert suffix == "grm.bin", f"Count mismatch: {first_path}"
		assert_floats_close(
			read_floats(first_path), read_floats(second_path), first_path, rel_tol, abs_tol
		)
	assert Path(f"{first}.grm.id").read_bytes() == Path(f"{second}.grm.id").read_bytes()
