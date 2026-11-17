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
