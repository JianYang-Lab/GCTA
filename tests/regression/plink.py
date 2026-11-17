import math
from pathlib import Path


def write_plink(prefix, genotypes):
	prefix = Path(prefix)
	prefix.parent.mkdir(parents=True, exist_ok=True)
	marker_count = len(genotypes)
	assert marker_count > 0
	sample_count = len(genotypes[0])
	assert sample_count > 0
	assert all(len(marker) == sample_count for marker in genotypes)
	assert all(value in (0, 1, 2, None) for marker in genotypes for value in marker)

	sample_ids = [("F1", f"I{index}") for index in range(sample_count)]
	prefix.with_suffix(".fam").write_text(
		"".join(f"{fid} {iid} 0 0 0 -9\n" for fid, iid in sample_ids), encoding="utf-8"
	)
	prefix.with_suffix(".bim").write_text(
		"".join(f"1 rs{index} 0 {1000 + index} A G\n" for index in range(marker_count)),
		encoding="utf-8",
	)
	prefix.with_suffix(".map").write_text(
		"".join(f"1 rs{index} 0 {1000 + index}\n" for index in range(marker_count)),
		encoding="utf-8",
	)

	allele_pairs = {0: ("G", "G"), 1: ("A", "G"), 2: ("A", "A"), None: ("0", "0")}
	with prefix.with_suffix(".ped").open("w", encoding="utf-8") as stream:
		for sample_index, (fid, iid) in enumerate(sample_ids):
			alleles = [
				allele
				for marker in genotypes
				for allele in allele_pairs[marker[sample_index]]
			]
			stream.write(f"{fid} {iid} 0 0 0 -9 " + " ".join(alleles) + "\n")

	code_by_genotype = {2: 0, None: 1, 1: 2, 0: 3}
	with prefix.with_suffix(".bed").open("wb") as stream:
		stream.write(b"\x6c\x1b\x01")
		for marker in genotypes:
			for offset in range(0, sample_count, 4):
				value = 0
				for within_byte, genotype in enumerate(marker[offset:offset + 4]):
					value |= code_by_genotype[genotype] << (2 * within_byte)
				stream.write(bytes((value,)))
	return sample_ids


def expected_grm(genotypes, algorithm):
	assert algorithm in (0, 1)
	assert all(value is not None for marker in genotypes for value in marker), (
		"The independent GRM oracle currently requires complete genotypes"
	)
	sample_count = len(genotypes[0])
	element_count = sample_count * (sample_count + 1) // 2
	sums = [0.0] * element_count
	variance_sum = 0.0
	valid_markers = 0

	for marker in genotypes:
		frequency = sum(marker) / (2.0 * sample_count)
		variance = 2.0 * frequency * (1.0 - frequency)
		if variance < 1e-50:
			continue
		valid_markers += 1
		variance_sum += variance
		if algorithm == 0:
			scale = math.sqrt(1.0 / variance)
			values = [(value - 2.0 * frequency) * scale for value in marker]
		else:
			values = [value - 2.0 * frequency for value in marker]
		for row in range(sample_count):
			offset = row * (row + 1) // 2
			for column in range(row + 1):
				sums[offset + column] += values[row] * values[column]

	assert valid_markers > 0
	if algorithm == 0:
		grm = [value / valid_markers for value in sums]
	else:
		grm = [value / variance_sum for value in sums]
	counts = [float(valid_markers)] * element_count
	return grm, counts


def synthetic_markers(sample_count=7, marker_count=131):
	markers = [
		[(sample * 3 + marker * 5 + sample * marker) % 3 for sample in range(sample_count)]
		for marker in range(marker_count)
	]
	markers[0] = [0] * sample_count
	markers[-1] = [2] * sample_count
	return markers
