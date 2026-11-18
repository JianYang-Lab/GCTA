from pathlib import Path
import struct

from plink import write_plink


def write_pgen(prefix, genotypes, chromosomes=None, dosage=False):
    prefix = Path(prefix)
    sample_ids = write_plink(prefix, genotypes, chromosomes)
    if chromosomes is None:
        chromosomes = [1] * len(genotypes)
    code_by_genotype = {0: 0, 1: 1, 2: 2, None: 3}
    records = bytearray()
    for marker in genotypes:
        for offset in range(0, len(sample_ids), 4):
            value = 0
            for within_byte, genotype in enumerate(marker[offset:offset + 4]):
                value |= code_by_genotype[genotype] << (2 * within_byte)
            records.append(value)
        if dosage:
            for genotype in marker:
                records.extend(struct.pack("<H", 0 if genotype is None else genotype * 16384))
    pgen = (
        b"\x6c\x1b"
        + bytes((3 if dosage else 2,))
        + struct.pack("<II", len(genotypes), len(sample_ids))
        + b"\x00"
        + records
    )
    prefix.with_suffix(".pgen").write_bytes(pgen)
    prefix.with_suffix(".pvar").write_text(
        "#CHROM\tID\tCM\tPOS\tALT\tREF\n"
        + "".join(
            f"{chromosomes[index]}\trs{index}\t0\t{1000 + index}\tA\tG\n"
            for index in range(len(genotypes))
        ),
        encoding="utf-8",
    )
    prefix.with_suffix(".psam").write_text(
        "#FID\tIID\tPAT\tMAT\tSEX\n"
        + "".join(f"{fid}\t{iid}\t0\t0\t0\n" for fid, iid in sample_ids),
        encoding="utf-8",
    )
    return sample_ids
