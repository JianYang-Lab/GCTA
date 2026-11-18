#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from gcta_runner import assert_floats_close, read_floats, run_gcta
from plink import synthetic_markers

MARKERS = synthetic_markers(sample_count=8, marker_count=4)


def expected_grm():
    sample_count = len(MARKERS[0])
    values = []
    for row in range(sample_count):
        for column in range(row + 1):
            total = 0.0
            for marker in MARKERS:
                mean = sum(marker) / sample_count
                variance = sum((value - mean) ** 2 for value in marker) / (sample_count - 1)
                scale = math.sqrt(1.0 / variance) if variance > 1e-50 else 0.0
                total += (marker[row] - mean) * (marker[column] - mean) * scale * scale
            values.append(total / len(MARKERS))
    return values


def write_fixture(work_dir, format_name):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    dose = work_dir / f"{format_name}.dose"
    info = work_dir / f"{format_name}.info"
    if format_name == "mach":
        dose.write_text(
            "".join(
                f"F1:I{sample} DOSE " + " ".join(str(marker[sample]) for marker in MARKERS) + "\n"
                for sample in range(len(MARKERS[0]))
            ),
            encoding="utf-8",
        )
        info.write_text(
            "SNP A1 A2 0 0 0 Rsq\n"
            + "".join(f"rs{index} A G 0 0 0 0.9\n" for index in range(len(MARKERS))),
            encoding="utf-8",
        )
    else:
        dose.write_text(
            "0 0 0 " + " ".join(f"I{sample}" for sample in range(len(MARKERS[0]))) + "\n"
            + "".join(
                f"rs{index} A G " + " ".join(str(value) for value in marker) + "\n"
                for index, marker in enumerate(MARKERS)
            ),
            encoding="utf-8",
        )
        info.write_text(
            "#CHROM SNP BP A1 A2 X X X 0 0 0 0 0\n"
            + "".join(
                f"1 rs{index} {1000 + index} A G x x x 0 0 0 0 0\n"
                for index in range(len(MARKERS))
            ),
            encoding="utf-8",
        )
    return dose, info


def main():
    parser = argparse.ArgumentParser(description="GCTA dosage input regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("mach", "beagle"))
    args = parser.parse_args()
    expected = expected_grm()
    for threads in (1, 4):
        work = Path(args.work) / f"threads_{threads}"
        dose, info = write_fixture(work, args.scenario)
        output = run_gcta(
            args.gcta,
            work,
            "grm",
            [f"--dosage-{args.scenario}", dose, info, "--make-grm"],
            threads=threads,
        )
        assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected, output)
        assert read_floats(output.with_suffix(".grm.N.bin")) == [4.0] * 36
        ids = output.with_suffix(".grm.id").read_text(encoding="utf-8").splitlines()
        if args.scenario == "mach":
            assert ids == [f"F1:I{index}\tF1:I{index}" for index in range(8)]
            info_lines = Path(info).read_text(encoding="utf-8").splitlines()
            fields = info_lines[1].split()
            fields[-1] = "0.5"
            info_lines[1] = " ".join(fields)
            filtered_info = Path(work) / "filtered.info"
            filtered_info.write_text("\n".join(info_lines) + "\n", encoding="utf-8")
            filtered = run_gcta(
                args.gcta,
                work,
                "filtered",
                ["--dosage-mach", dose, filtered_info, "--imput-rsq", 0.8, "--make-grm"],
                threads=threads,
            )
            assert "there are 3 SNPs" in filtered.with_suffix(".log").read_text(encoding="utf-8")
        else:
            assert ids == [f"I{index}\tI{index}" for index in range(8)]
    print(f"Passed {args.scenario} dosage GRM, ID, and 1/4-thread regression checks")


if __name__ == "__main__":
    main()
