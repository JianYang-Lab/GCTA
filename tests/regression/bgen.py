import sqlite3
import struct
from pathlib import Path
import zlib


def write_bgen(prefix, genotypes, chromosomes=None, sample_indices=None):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    marker_count = len(genotypes)
    sample_count = len(genotypes[0])
    if chromosomes is None:
        chromosomes = [1] * marker_count
    assert all(len(marker) == sample_count for marker in genotypes)
    assert len(chromosomes) == marker_count

    variants = []
    bgen_path = prefix.with_suffix(".bgen")
    with bgen_path.open("wb") as stream:
        stream.write(struct.pack("<I", 20))
        stream.write(struct.pack("<III", 20, marker_count, sample_count))
        stream.write(b"bgen")
        stream.write(struct.pack("<I", 1 | (2 << 2)))
        for marker_index, marker in enumerate(genotypes):
            variant_start = stream.tell()
            rsid = f"rs{marker_index}".encode()
            chromosome = str(chromosomes[marker_index]).encode()
            stream.write(struct.pack("<H", 0))
            stream.write(struct.pack("<H", len(rsid)))
            stream.write(rsid)
            stream.write(struct.pack("<H", len(chromosome)))
            stream.write(chromosome)
            stream.write(struct.pack("<I H", 1000 + marker_index, 2))
            for allele in (b"A", b"G"):
                stream.write(struct.pack("<I", len(allele)))
                stream.write(allele)

            probabilities = bytearray()
            ploidy = bytearray()
            for genotype in marker:
                if genotype is None:
                    ploidy.append(0x80 | 2)
                    probabilities.extend((0, 0))
                else:
                    ploidy.append(2)
                    probabilities.extend({0: (0, 0), 1: (0, 255), 2: (255, 0)}[genotype])
            payload = struct.pack("<I H B B", sample_count, 2, 2, 2)
            payload += ploidy + b"\x00\x08" + probabilities
            compressed = zlib.compress(payload)
            stream.write(struct.pack("<II", len(compressed) + 4, len(payload)))
            stream.write(compressed)
            variants.append((
                chromosome.decode(), 1000 + marker_index, rsid.decode(), "A", "G",
                variant_start, stream.tell() - variant_start, 2,
            ))

    if sample_indices is None:
        sample_indices = range(sample_count)
    assert len(sample_indices) == sample_count
    with prefix.with_suffix(".sample").open("w", encoding="utf-8") as stream:
        stream.write("ID_1 ID_2 missing sex\n0 0 0 D\n")
        for sample in sample_indices:
            stream.write(f"F1 I{sample} 0 0\n")

    index_path = Path(f"{bgen_path}.bgi")
    if index_path.exists():
        index_path.unlink()
    connection = sqlite3.connect(index_path)
    with connection:
        connection.execute(
            "CREATE TABLE Metadata (filename TEXT, file_size INTEGER, reserved INTEGER, first_1000_bytes BLOB)"
        )
        connection.execute(
            "INSERT INTO Metadata VALUES (?, ?, 0, ?)",
            (bgen_path.name, bgen_path.stat().st_size, bgen_path.read_bytes()[:1000]),
        )
        connection.execute(
            "CREATE TABLE Variant (chromosome TEXT, position INTEGER, rsid TEXT, allele1 TEXT, allele2 TEXT, file_start_position INTEGER, size_in_bytes INTEGER, number_of_alleles INTEGER)"
        )
        connection.executemany("INSERT INTO Variant VALUES (?, ?, ?, ?, ?, ?, ?, ?)", variants)
    connection.close()
