#!/usr/bin/env python3

import argparse
from pathlib import Path

from gcta_runner import run_gcta
from test_errors import expect_error

GENOTYPES = [
    [("A", "G", 0.95), ("C", "T", 0.69), ("A", "A", 0.70)],
    [("G", "G", 0.80), ("C", "C", 0.90), ("G", "G", 0.60)],
]


def write_raw(path, version, sample, mismatch=False):
    version = "1.1.9" if version == "old" else "2.0.0"
    header = ["[Header]"] * 11
    header[1] = f"[Header] Version {version}"
    header[4] = "Total SNPs 3"
    rows = []
    for index, (a1, a2, gencall) in enumerate(GENOTYPES[sample]):
        snp = f"rs{index}"
        if version == "1.1.9":
            rows.append(f"{snp} I{sample} {a1} {a2} 0 0 0 0 {gencall}\n")
        else:
            rows.append(f"{snp} I{sample} F{sample} {gencall} 0 0 {a1} {a2}\n")
    if mismatch:
        rows[:2] = [rows[1], rows[0]]
    Path(path).write_text("\n".join(header) + "\n" + "".join(rows), encoding="utf-8")


def make_fixture(work_dir, version, mismatch=False):
    work_dir = Path(work_dir)
    raw_dir = work_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_files = []
    for sample in range(2):
        raw_file = raw_dir / f"sample{sample}.txt"
        write_raw(raw_file, version, sample, mismatch and sample == 0)
        raw_files.append(raw_file)
    file_list = work_dir / "files.txt"
    file_list.write_text("".join(f"{path}\n" for path in raw_files), encoding="utf-8")
    summary = work_dir / "summary.txt"
    summary.write_text(
        "Index SNP CHR BP\n0 rs0 1 1000\n1 rs1 1 2000\n2 rs2 1 3000\n",
        encoding="utf-8",
    )
    return file_list, summary


def expected_ped(version):
    lines = []
    for sample, markers in enumerate(GENOTYPES):
        fid = f"I{sample}" if version == "old" else f"F{sample}"
        alleles = []
        for a1, a2, gencall in markers:
            if gencall < 0.7:
                alleles.extend(("0", "0"))
            else:
                alleles.extend((a1, a2))
        lines.append(f"{fid} I{sample} -9 -9 -9 -9 " + " ".join(alleles) + " \n")
    return "".join(lines)


def run_conversion(gcta, work_dir, version, threads):
    file_list, summary = make_fixture(work_dir, version)
    output = run_gcta(
        gcta,
        work_dir,
        "raw",
        ["--raw-files", file_list, "--raw-summary", summary],
        threads=threads,
    )
    assert output.with_suffix(".ped").read_text(encoding="utf-8") == expected_ped(version)
    assert output.with_suffix(".map").read_text(encoding="utf-8") == (
        "1\trs0\t0\t1000\n1\trs1\t0\t2000\n1\trs2\t0\t3000\n"
    )
    return output


def main():
    parser = argparse.ArgumentParser(description="GCTA raw genotype conversion regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=("modern", "old", "mismatch"))
    args = parser.parse_args()
    if args.scenario == "mismatch":
        work = Path(args.work)
        file_list, summary = make_fixture(work, "modern", mismatch=True)
        output = work / "mismatch"
        expect_error(
            args.gcta,
            work,
            "mismatch",
            ["--raw-files", file_list, "--raw-summary", summary, "--out", output],
        )
        console = (work / "mismatch.console").read_text(encoding="utf-8")
        assert "Has the order of the SNPs been changed?" in console
    else:
        single = run_conversion(args.gcta, Path(args.work) / "single", args.scenario, 1)
        multi = run_conversion(args.gcta, Path(args.work) / "multi", args.scenario, 4)
        assert single.with_suffix(".ped").read_bytes() == multi.with_suffix(".ped").read_bytes()
        assert single.with_suffix(".map").read_bytes() == multi.with_suffix(".map").read_bytes()
    print(f"Passed raw genotype regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
