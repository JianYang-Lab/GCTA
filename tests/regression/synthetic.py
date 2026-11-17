import math
from pathlib import Path
import random


def family_markers(sample_count=128, marker_count=2049):
	assert sample_count >= 4 and sample_count % 2 == 0
	assert marker_count >= 3
	founder_count = sample_count // 2
	rng = random.Random(20261117)
	markers = []
	chromosomes = []

	for marker_index in range(marker_count):
		frequency = 0.02 + (marker_index % 480) / 1000.0
		marker = [rng.randrange(2) + rng.randrange(2) for _ in range(founder_count)]
		for offspring in range(founder_count):
			family = offspring // 2
			parents = (2 * family, 2 * family + 1)
			genotype = 0
			for parent in parents:
				parent_genotype = marker[parent]
				if parent_genotype == 2:
					allele = 1
				elif parent_genotype == 0:
					allele = 0
				else:
					allele = rng.randrange(2)
				genotype += allele
			marker.append(genotype)

		if marker_index == 0:
			marker = [0] * sample_count
		elif marker_index == marker_count - 1:
			marker = [2] * sample_count
		else:
			for sample in range(sample_count):
				if (marker_index * 13 + sample * 17) % 997 == 0:
					marker[sample] = None
		markers.append(marker)
		chromosomes.append(1 + marker_index * 3 // marker_count)
	return markers, chromosomes


def write_family_phenotypes(prefix, markers, case_fraction=0.5):
	prefix = Path(prefix)
	sample_count = len(markers[0])
	marker_count = len(markers)
	rng = random.Random(20261118)
	causal_indices = [16 + index * (marker_count - 32) // 32 for index in range(32)]
	genetic = [0.0] * sample_count
	for effect_index, marker_index in enumerate(causal_indices):
		marker = markers[marker_index]
		observed = [value for value in marker if value is not None]
		frequency = sum(observed) / (2.0 * len(observed))
		variance = 2.0 * frequency * (1.0 - frequency)
		if variance < 1e-50:
			continue
		scale = math.sqrt(1.0 / variance)
		effect = (0.5 + (effect_index % 5) / 5.0) * (1 if effect_index % 2 == 0 else -1)
		for sample in range(sample_count):
			if marker[sample] is not None:
				genetic[sample] += effect * (marker[sample] - 2.0 * frequency) * scale

	mean = sum(genetic) / sample_count
	variance = sum((value - mean) ** 2 for value in genetic) / sample_count
	assert variance > 0
	genetic = [(value - mean) / math.sqrt(variance) for value in genetic]
	quantitative_covariates = [
		(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)) for _ in range(sample_count)
	]
	quantitative = []
	for sample in range(sample_count):
		qcovar1, qcovar2 = quantitative_covariates[sample]
		value = (
			math.sqrt(0.7) * genetic[sample]
			+ math.sqrt(0.3) * rng.gauss(0.0, 1.0)
			+ 0.2 * qcovar1
			- 0.1 * qcovar2
		)
		quantitative.append(value)
	assert 0.0 < case_fraction < 1.0
	case_count = max(1, min(sample_count - 1, round(sample_count * case_fraction)))
	ordered = sorted(quantitative)
	threshold_index = sample_count - case_count
	threshold = (ordered[threshold_index - 1] + ordered[threshold_index]) / 2.0
	binary = [1 if value > threshold else 0 for value in quantitative]

	prefix.with_suffix(".phen").write_text(
		"".join(
			f"F1\tI{sample}\t{quantitative[sample]:.17g}\t{binary[sample]}\n"
			for sample in range(sample_count)
		),
		encoding="utf-8",
	)
	prefix.with_suffix(".qcovar").write_text(
		"".join(
			f"F1\tI{sample}\t{qcovar1:.17g}\t{qcovar2:.17g}\n"
			for sample, (qcovar1, qcovar2) in enumerate(quantitative_covariates)
		),
		encoding="utf-8",
	)
	return prefix.with_suffix(".phen"), prefix.with_suffix(".qcovar")
