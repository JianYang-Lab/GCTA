#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, run_gcta
from plink import write_plink

MARKERS = [
    [0, 1, 2, 0, 1, 2],
    [0, 0, 1, 1, 2, 2],
    [2, 2, 2, 0, 1, 0],
    [None, 1, None, 0, 2, 1],
]
SEXES = [1, 2, 1, 2, 0, 1]


def make_fixture(work_dir, chromosome):
    prefix = Path(work_dir) / "input" / "x"
    write_plink(prefix, MARKERS, chromosomes=[chromosome] * len(MARKERS))
    lines = prefix.with_suffix(".fam").read_text(encoding="utf-8").splitlines()
    prefix.with_suffix(".fam").write_text(
        "".join(
            f"{line.split()[0]}\t{line.split()[1]}\t0\t0\t{sex}\t-9\n"
            for line, sex in zip(lines, SEXES)
        ),
        encoding="utf-8",
    )
    return prefix


def marker_statistics(marker, x_mode):
    valid_n = 0
    male_valid_n = 0
    dosage = 0.0
    male_hom_alt_n = 0
    male_het_n = 0
    for value, sex in zip(marker, SEXES):
        if value is None:
            continue
        valid_n += 1
        dosage += value
        if sex == 1:
            male_valid_n += 1
            if value == 2:
                male_hom_alt_n += 1
            elif value == 1:
                male_het_n += 1
    if x_mode:
        valid_alleles = 2 * valid_n - male_valid_n
        af = (dosage - male_hom_alt_n - 0.5 * male_het_n) / valid_alleles
        mean = 2.0 * af
    else:
        af = dosage / (2 * valid_n)
        mean = 2.0 * af
    return af, mean, mean * (1.0 - af)


def expected_grm(x_mode, dominance):
    sample_count = len(MARKERS[0])
    transformed = []
    valid_counts = []
    for marker in MARKERS:
        af, mean, sd = marker_statistics(marker, x_mode)
        if sd < 1e-50:
            transformed.append([0.0] * sample_count)
            valid_counts.append(False)
            continue
        values = []
        for index, (value, sex) in enumerate(zip(marker, SEXES)):
            if dominance:
                center = 0.5 * mean * mean
                scale = 1.0 / sd
                if value is None:
                    item = (center - center) * scale
                else:
                    allele_value = {0: 0.0, 1: mean, 2: 2.0 * mean - 2.0}[value]
                    item = (allele_value - center) * scale
            else:
                center = mean if x_mode else 2.0 * af
                scale = math.sqrt(1.0 / sd)
                item = mean - center if value is None else (value - center) * scale
                if x_mode and sex == 1:
                    item *= math.sqrt(0.5)
            values.append(item)
        transformed.append(values)
        valid_counts.append(True)

    expected = []
    counts = []
    for row in range(sample_count):
        for column in range(row + 1):
            total = 0.0
            count = 0
            for index, marker in enumerate(MARKERS):
                if not valid_counts[index]:
                    continue
                if MARKERS[index][row] is None or MARKERS[index][column] is None:
                    continue
                total += transformed[index][row] * transformed[index][column]
                count += 1
            expected.append(total / count if count else 0.0)
            counts.append(float(count))
    return expected, counts


def main():
    parser = argparse.ArgumentParser(description="GCTA X and dominance GRM regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("x", "dominance"))
    args = parser.parse_args()
    expected, expected_counts = expected_grm(args.scenario == "x", args.scenario == "dominance")
    outputs = []
    for threads in (1, 4):
        work = Path(args.work) / f"threads_{threads}"
        prefix = make_fixture(work, 23 if args.scenario == "x" else 1)
        option = "--make-grm-xchr" if args.scenario == "x" else "--make-grm-d"
        output = run_gcta(args.gcta, work, "grm", ["--bfile", prefix, option], threads=threads)
        suffix = ".d.grm.bin" if args.scenario == "dominance" else ".grm.bin"
        assert_floats_close(read_floats(output.with_suffix(suffix)), expected, output)
        assert_floats_close(read_floats(output.with_suffix(suffix.replace(".grm.bin", ".grm.N.bin"))), expected_counts, output)
        outputs.append(output)
    if args.scenario == "dominance":
        for suffix in (".d.grm.bin", ".d.grm.N.bin", ".d.grm.id"):
            assert Path(f"{outputs[0]}{suffix}").read_bytes() == Path(f"{outputs[1]}{suffix}").read_bytes()
    else:
        assert_grm_equal(*outputs)
    print(f"Passed GRM scenario: {args.scenario}")


if __name__ == "__main__":
    main()
