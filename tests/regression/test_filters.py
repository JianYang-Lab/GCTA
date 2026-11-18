#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_floats_close, assert_grm_equal, read_floats, read_table, run_gcta
from plink import expected_grm, write_plink

MARKERS = [
    [0, 0, 0, 0, 2, 2, 2, 2],
    [0, 0, 0, 0, 0, 0, 2, 2],
    [0, 0, 0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1, 2, 2],
    [2, 2, 2, 2, 0, 0, 0, 0],
    [None, None, 0, 0, 0, 0, 1, 2],
    [0, 1, 2, 0, 1, 2, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0],
]
CHROMOSOMES = [1, 1, 1, 2, 2, 3, 3, 1]
SCENARIOS = (
    "maf",
    "max_maf",
    "chr",
    "chr_range",
    "autosome",
    "autosome_num",
    "remove",
    "exclude",
    "combined",
)
REMOVED_SAMPLES = {0, 7}


def marker_frequency(marker):
    observed = [value for value in marker if value is not None]
    frequency = sum(observed) / (2.0 * len(observed))
    return frequency, min(frequency, 1.0 - frequency), 2 * len(observed)


def scenario_configuration(scenario, work_dir):
    arguments = []
    marker_indices = set(range(len(MARKERS)))
    sample_indices = list(range(len(MARKERS[0])))
    minimum_maf = 0.0
    maximum_maf = 0.5
    filter_maf = False
    if scenario == "maf":
        arguments = ["--maf", 0.2]
        minimum_maf = 0.2
        filter_maf = True
    elif scenario == "max_maf":
        arguments = ["--max-maf", 0.3]
        maximum_maf = 0.3
        filter_maf = True
    elif scenario == "chr":
        arguments = ["--chr", 2]
        marker_indices = {index for index in marker_indices if CHROMOSOMES[index] == 2}
    elif scenario == "chr_range":
        arguments = ["--chr", 1, 2]
        marker_indices = {index for index in marker_indices if CHROMOSOMES[index] in (1, 2)}
    elif scenario == "autosome":
        arguments = ["--autosome"]
        marker_indices = {index for index in marker_indices if CHROMOSOMES[index] <= 22}
    elif scenario == "autosome_num":
        arguments = ["--autosome-num", 23]
    elif scenario == "remove":
        remove = work_dir / "remove.txt"
        remove.write_text("".join(f"F1 I{index}\n" for index in sorted(REMOVED_SAMPLES)), encoding="utf-8")
        arguments = ["--remove", remove]
        sample_indices = [index for index in sample_indices if index not in REMOVED_SAMPLES]
    elif scenario == "exclude":
        exclude = work_dir / "exclude.txt"
        exclude.write_text("rs1\nrs5\n", encoding="utf-8")
        arguments = ["--exclude", exclude]
        marker_indices -= {1, 5}
    elif scenario == "combined":
        exclude = work_dir / "combined_exclude.txt"
        exclude.write_text("rs1\n", encoding="utf-8")
        remove = work_dir / "combined_remove.txt"
        remove.write_text("".join(f"F1 I{index}\n" for index in sorted(REMOVED_SAMPLES)), encoding="utf-8")
        arguments = ["--chr", 1, 2, "--exclude", exclude, "--remove", remove, "--maf", 0.2]
        marker_indices = {index for index in marker_indices if CHROMOSOMES[index] in (1, 2)}
        marker_indices -= {1}
        sample_indices = [index for index in sample_indices if index not in REMOVED_SAMPLES]
        minimum_maf = 0.2
        filter_maf = True

    markers = [
        [MARKERS[index][sample] for sample in sample_indices]
        for index in sorted(marker_indices)
    ]
    if filter_maf:
        selected = [
            (index, marker)
            for index, marker in zip(sorted(marker_indices), markers)
            if minimum_maf <= marker_frequency(marker)[1] <= maximum_maf
        ]
        marker_indices = {index for index, _ in selected}
        markers = [marker for _, marker in selected]
    return arguments, sorted(marker_indices), sample_indices, markers


def check_frequency(path, marker_indices, sample_indices):
    rows = read_table(path)
    assert rows[0] == ["CHR", "SNP", "POS", "A1", "A2", "AF", "NCHROBS"]
    assert len(rows) == len(marker_indices) + 1
    for fields, index in zip(rows[1:], marker_indices):
        marker = [MARKERS[index][sample] for sample in sample_indices]
        frequency, _, chromosome_count = marker_frequency(marker)
        assert fields[:5] == [str(CHROMOSOMES[index]), f"rs{index}", str(1000 + index), "A", "G"]
        assert math.isclose(float(fields[5]), frequency, rel_tol=0.0, abs_tol=1e-6)
        assert fields[6] == str(chromosome_count)


def check_scenario(gcta, work_dir, scenario, threads):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    prefix = Path(work_dir.parent) / "input" / "filters"
    arguments, marker_indices, sample_indices, markers = scenario_configuration(scenario, work_dir)
    frequency = run_gcta(
        gcta,
        work_dir,
        "frequency",
        ["--bfile", prefix, *arguments, "--freq"],
        threads=threads,
    )
    check_frequency(frequency.with_suffix(".frq"), marker_indices, sample_indices)

    grm = run_gcta(
        gcta,
        work_dir,
        "grm",
        ["--bfile", prefix, *arguments, "--make-grm"],
        threads=threads,
    )
    expected_values, expected_counts = expected_grm(
        [
            marker
            for index, marker in zip(marker_indices, markers)
            if CHROMOSOMES[index] <= (23 if scenario == "autosome_num" else 22)
        ],
        0,
    )
    assert_floats_close(read_floats(grm.with_suffix(".grm.bin")), expected_values, grm)
    assert read_floats(grm.with_suffix(".grm.N.bin")) == expected_counts
    expected_ids = "".join(f"F1\tI{index}\n" for index in sample_indices)
    assert grm.with_suffix(".grm.id").read_text(encoding="utf-8") == expected_ids
    return frequency.with_suffix(".frq"), grm


def main():
    parser = argparse.ArgumentParser(description="GCTA genotype filtering regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    args = parser.parse_args()
    global CHROMOSOMES
    if args.scenario in ("autosome", "autosome_num"):
        CHROMOSOMES = CHROMOSOMES.copy()
        CHROMOSOMES[6] = 23
    work = Path(args.work)
    write_plink(work / "input" / "filters", MARKERS, CHROMOSOMES)
    fam = work / "input" / "filters.fam"
    fam.write_text(
        "".join(
            f"{fields[0]} {fields[1]} 0 0 2 -9\n"
            for fields in (line.split() for line in fam.read_text(encoding="utf-8").splitlines())
        ),
        encoding="utf-8",
    )

    single_frequency, single_grm = check_scenario(args.gcta, work / "single", args.scenario, 1)
    multi_frequency, multi_grm = check_scenario(args.gcta, work / "multi", args.scenario, 4)
    assert single_frequency.read_bytes() == multi_frequency.read_bytes()
    assert_grm_equal(single_grm, multi_grm)
    print(f"Passed genotype filtering scenario: {args.scenario}")


if __name__ == "__main__":
    main()
