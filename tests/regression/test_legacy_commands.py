#!/usr/bin/env python3

import argparse
import gzip
import math
from pathlib import Path
import random
import struct

from gcta_runner import (
    assert_floats_close,
    assert_grm_equal,
    assert_tables_close,
    read_floats,
    read_table,
    run_gcta,
)
from plink import synthetic_markers, write_plink
from synthetic import family_markers, write_family_phenotypes
from test_grm_x import MARKERS as X_MARKERS
from test_grm_x import SEXES as X_SEXES
from test_grm_x import expected_grm as expected_special_grm
from test_ld_legacy import MARKERS as LD_MARKERS
from test_legacy import expected_legacy_grm
from test_raw import GENOTYPES as RAW_GENOTYPES
from test_raw import make_fixture as make_raw_fixture


SCENARIOS = (
    "raw-files",
    "hapmap-genet-dst",
    "paa",
    "ibc",
    "make-grm-formats",
    "make-grm-d-formats",
    "make-grm-xchr-formats",
    "make-grm-inbred-formats",
    "recode",
    "ld",
    "ld-score-multi",
    "ld-score-region",
    "ld-max-rsq",
    "blup-snp",
    "cojo-select",
    "cojo-cond",
    "cojo-sblup",
    "fastbat-set",
    "fastbat-seg",
    "mbat-set",
    "gwas-adj-pc",
    "pc-loading",
    "project-loading",
    "hereg-bivar",
    "reml-variants",
    "grm-align",
    "make-bk",
)


def path_with(prefix, suffix):
    return Path(f"{prefix}{suffix}")


def gzip_table(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line.split() for line in stream.read().splitlines()]


def assert_gzip_tables_close(actual, expected):
    actual_rows = gzip_table(actual)
    expected_rows = gzip_table(expected)
    assert len(actual_rows) == len(expected_rows)
    for row_index, (actual_row, expected_row) in enumerate(zip(actual_rows, expected_rows)):
        assert len(actual_row) == len(expected_row), (actual, row_index)
        for column_index, (observed, reference) in enumerate(zip(actual_row, expected_row)):
            if observed == reference:
                continue
            try:
                observed_number = float(observed)
                reference_number = float(reference)
            except ValueError as error:
                raise AssertionError(
                    f"Text mismatch at {actual}:{row_index + 1}:{column_index + 1}: "
                    f"{observed} != {reference}"
                ) from error
            assert math.isclose(observed_number, reference_number, rel_tol=3e-5, abs_tol=2e-6), (
                f"Value mismatch at {actual}:{row_index + 1}:{column_index + 1}: "
                f"{observed} != {reference}"
            )


def run_legacy(gcta, work, name, arguments, threads, terminal_log):
    output = run_gcta(gcta, work, name, arguments, threads=threads, timeout=120)
    log = output.with_suffix(".log").read_text(encoding="utf-8")
    assert terminal_log in log, f"Legacy terminal path not reached for {name}: {terminal_log}\n{log}"
    assert "Accepted options:" in log
    return output


def assert_close(observed, expected, context="", rel_tol=3e-5, abs_tol=2e-6):
    assert math.isfinite(observed), f"Non-finite value for {context}: {observed}"
    assert math.isclose(observed, expected, rel_tol=rel_tol, abs_tol=abs_tol), (
        f"Value mismatch for {context}: {observed} != {expected}"
    )


def make_basic_fixture(work, sample_count=8, marker_count=12):
    markers = synthetic_markers(sample_count, marker_count)
    markers[1][0] = None
    prefix = Path(work) / "input" / "synthetic"
    write_plink(prefix, markers)
    return prefix, markers


def make_polymorphic_fixture(work):
    rng = random.Random(4401)
    sample_count = 8
    markers = []
    for _ in range(10):
        while True:
            marker = [rng.randrange(3) for _ in range(sample_count)]
            if len(set(marker)) == 3:
                break
        markers.append(marker)
    markers[1][0] = None
    prefix = Path(work) / "input" / "polymorphic"
    write_plink(prefix, markers)
    return prefix, markers


def compare_suffixes(single, multi, suffixes):
    for suffix in suffixes:
        assert_tables_close(path_with(single, suffix), path_with(multi, suffix))


def scenario_raw_files(gcta, work):
    for version in ("modern", "old"):
        outputs = []
        for threads in (1, 4):
            case = Path(work) / version / f"threads_{threads}"
            file_list, summary = make_raw_fixture(case, version)
            output = run_legacy(
                gcta,
                case,
                "raw",
                ["--raw-files", file_list, "--raw-summary", summary, "--gencall", 0.8],
                threads,
                "Reading the raw genotype files and saving the genotype data in PLINK PED format",
            )
            expected = []
            for sample, markers in enumerate(RAW_GENOTYPES):
                fid = f"I{sample}" if version == "old" else f"F{sample}"
                alleles = []
                for a1, a2, gencall in markers:
                    alleles.extend(("0", "0") if gencall < 0.8 else (a1, a2))
                expected.append(f"{fid} I{sample} -9 -9 -9 -9 " + " ".join(alleles) + " \n")
            assert output.with_suffix(".ped").read_text(encoding="utf-8") == "".join(expected)
            assert output.with_suffix(".map").read_text(encoding="utf-8") == (
                "1\trs0\t0\t1000\n1\trs1\t0\t2000\n1\trs2\t0\t3000\n"
            )
            outputs.append(output)
        for suffix in (".ped", ".map"):
            assert path_with(outputs[0], suffix).read_bytes() == path_with(outputs[1], suffix).read_bytes()


def scenario_hapmap(gcta, work):
    prefix, markers = make_basic_fixture(work)
    map_prefix = Path(work) / "input" / "genetic_map_"
    for chromosome in range(1, 23):
        Path(f"{map_prefix}{chromosome}_CEU_b36.txt").write_text(
            "Position Rate(cM/Mb) Map(cM)\n0 1 0\n1000000 1 1000\n", encoding="utf-8"
        )
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "distance",
            ["--bfile", prefix, "--hapmap-genet-dst", map_prefix],
            threads,
            "Genetic distances have been created",
        )
        rows = read_table(path_with(output, ".genetdst"))
        assert len(rows) == len(markers)
        for index, row in enumerate(rows):
            assert row[:2] == ["1", f"rs{index}"]
            assert_close(float(row[2]), 0.0 if index == 0 else 1e-6, f"genetic distance {index}")
            assert row[3:] == [str(1000 + index), "A", "G"]
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".genetdst",))


def expected_paa(markers, ancestral):
    sample_count = len(markers[0])
    expected = []
    for sample in range(sample_count):
        counts = [0.0] * 6
        nomiss = 0
        for marker, ancestor in zip(markers, ancestral):
            observed = [value for value in marker if value is not None]
            mu = sum(observed) / len(observed)
            value = marker[sample]
            if value is None or ancestor == ".":
                continue
            nomiss += 1
            if value == 0:
                if ancestor == "A":
                    category = 2 if mu > 1.0 else 3
                else:
                    category = 0 if mu > 1.0 else 1
            elif value == 2:
                if ancestor == "A":
                    category = 1 if mu > 1.0 else 0
                else:
                    category = 3 if mu > 1.0 else 2
            else:
                if ancestor == "A":
                    category = 5 if mu > 1.0 else 4
                else:
                    category = 4 if mu > 1.0 else 5
            counts[category] += 1.0
        expected.append([float(nomiss)] + [value / nomiss for value in counts])
    return expected


def scenario_paa(gcta, work):
    prefix, markers = make_basic_fixture(work)
    ancestral = ["A" if index % 2 else "G" for index in range(len(markers))]
    ancestral_file = Path(work) / "input" / "ancestral.txt"
    ancestral_file.write_text(
        "".join(f"rs{index} {allele}\n" for index, allele in enumerate(ancestral)), encoding="utf-8"
    )
    expected = expected_paa(markers, ancestral)
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "paa",
            ["--bfile", prefix, "--paa", ancestral_file],
            threads,
            "Calculating proportion of ancestral alleles",
        )
        rows = read_table(output.with_suffix(".paa"))
        assert rows[0] == [
            "FID", "IID", "NOMISS", "HOM_AA_RARE", "HOM_AA_COMM", "HOM_DA_RARE",
            "HOM_DA_COMM", "HET_AA_RARE", "HET_AA_COMM",
        ]
        assert len(rows) == len(markers[0]) + 1
        for sample, (row, reference) in enumerate(zip(rows[1:], expected)):
            assert row[:2] == ["F1", f"I{sample}"]
            for column, value in enumerate(reference):
                assert_close(float(row[column + 2]), value, f"PAA sample {sample} column {column}")
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".paa",))


def expected_ibc(markers, all_estimators):
    sample_count = len(markers[0])
    rows = []
    for sample in range(sample_count):
        sums = [0.0] * 11
        nomiss = 0.0
        sum_w = 0.0
        sum_h = 0.0
        for marker in markers:
            observed = [value for value in marker if value is not None]
            mu = sum(observed) / len(observed)
            p = 0.5 * mu
            h = 2.0 * p * (1.0 - p)
            h_i = 0.0 if abs(h) < 1e-50 else 1.0 / h
            w = h / (1.0 - h)
            p_q = p / (1.0 - p) if p != 1.0 else math.inf
            q_p = 0.0 if abs(p_q) < 1e-50 or math.isinf(p_q) else 1.0 / p_q
            value = marker[sample]
            if value is None:
                continue
            x = float(value)
            base = (x - mu) ** 2
            f1 = base * h_i
            f2 = (h - x * (2.0 - x)) * h_i
            f3 = (x * (x - 1.0 - mu) + mu * p) * h_i
            sums[0] += f1
            sums[2] += f2
            sums[4] += f3
            if all_estimators:
                sums[1] += f1 * w
                sums[3] += f2 * w
                sums[5] += base
                sums[6] += x * (x - 1.0 - mu) + mu * p
                sums[7] += h - x * (2.0 - x)
                if x < 0.1:
                    sums[8] += p_q
                    if p > 0.5:
                        sums[9] += 1.0
                    else:
                        sums[10] += 1.0
                elif x > 1.9:
                    sums[8] += q_p
                    if p < 0.5:
                        sums[9] += 1.0
                    else:
                        sums[10] += 1.0
                sum_w += w
                sum_h += h
            nomiss += 1.0
        simple = [nomiss, sums[0] / nomiss - 1.0, sums[2] / nomiss, sums[4] / nomiss]
        if not all_estimators:
            rows.append(simple)
            continue
        rows.append([
            nomiss,
            sums[9] / nomiss,
            sums[10] / nomiss,
            simple[1],
            sums[1] / sum_w - 1.0,
            simple[2],
            sums[3] / sum_w,
            simple[3],
            sums[5] / sum_h - 1.0,
            sums[6] / sum_h,
            sums[7] / sum_h,
            sums[8] / nomiss,
        ])
    return rows


def scenario_ibc(gcta, work):
    prefix, markers = make_basic_fixture(work)
    for option in ("--ibc", "--ibc-all"):
        all_estimators = option == "--ibc-all"
        expected = expected_ibc(markers, all_estimators)
        outputs = []
        for threads in (1, 4):
            output = run_legacy(
                gcta,
                Path(work) / option[2:] / f"threads_{threads}",
                "ibc",
                ["--bfile", prefix, option],
                threads,
                "Calculating the inbreeding coefficients",
            )
            rows = read_table(output.with_suffix(".ibc"))
            assert len(rows) == len(markers[0]) + 1
            assert len(rows[0]) == (14 if all_estimators else 6)
            for sample, (row, reference) in enumerate(zip(rows[1:], expected)):
                assert row[:2] == ["F1", f"I{sample}"]
                for column, value in enumerate(reference):
                    assert_close(float(row[column + 2]), value, f"IBC sample {sample} column {column}")
            outputs.append(output)
        compare_suffixes(outputs[0], outputs[1], (".ibc",))


def assert_gz_grm(output, stem, expected_values, expected_counts):
    rows = gzip_table(path_with(output, f"{stem}.grm.gz"))
    assert len(rows) == len(expected_values)
    sample_count = len(read_table(path_with(output, f"{stem}.grm.id")))
    expected_indices = [(row, column) for row in range(sample_count) for column in range(row + 1)]
    for index, (row, value, count) in enumerate(zip(rows, expected_values, expected_counts)):
        assert [int(row[0]) - 1, int(row[1]) - 1] == list(expected_indices[index])
        assert_close(float(row[2]), count, f"GRM count {index}")
        assert_close(float(row[3]), value, f"GRM value {index}")


def scenario_make_grm_formats(gcta, work):
    prefix, markers = make_basic_fixture(work)
    expected, counts = expected_legacy_grm(markers, 0)
    outputs = {}
    for option in ("--make-grm-bin", "--make-grm-gz"):
        for threads in (1, 4):
            output = run_legacy(
                gcta,
                Path(work) / option[2:] / f"threads_{threads}",
                "grm",
                ["--bfile", prefix, option],
                threads,
                "Calculating the genetic relationship matrix",
            )
            if option.endswith("bin"):
                assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected, output)
                assert_floats_close(read_floats(output.with_suffix(".grm.N.bin")), counts, output)
            else:
                assert_gz_grm(output, "", expected, counts)
            outputs[(option, threads)] = output
        if option.endswith("bin"):
            assert_grm_equal(outputs[(option, 1)], outputs[(option, 4)])
        else:
            assert_gzip_tables_close(path_with(outputs[(option, 1)], ".grm.gz"), path_with(outputs[(option, 4)], ".grm.gz"))


def make_x_fixture(work, chromosome):
    prefix = Path(work) / "input" / "x"
    write_plink(prefix, X_MARKERS, chromosomes=[chromosome] * len(X_MARKERS))
    lines = prefix.with_suffix(".fam").read_text(encoding="utf-8").splitlines()
    prefix.with_suffix(".fam").write_text(
        "".join(
            f"{line.split()[0]}\t{line.split()[1]}\t0\t0\t{sex}\t-9\n"
            for line, sex in zip(lines, X_SEXES)
        ),
        encoding="utf-8",
    )
    return prefix


def scenario_special_grm(gcta, work, kind):
    x_mode = kind == "xchr"
    prefix = make_x_fixture(work, 23 if x_mode else 1)
    expected, counts = expected_special_grm(x_mode, not x_mode)
    if kind == "d":
        options = ("--make-grm-d-v1", "--make-grm-d-bin", "--make-grm-d-gz")
        stem = ".d"
    else:
        options = ("--make-grm-xchr-bin", "--make-grm-xchr-gz")
        stem = ""
    outputs = {}
    sex_file = None
    if x_mode:
        sex_file = Path(work) / "input" / "sex.txt"
        sex_file.write_text(
            "".join(
                f"F1 I{index} {2 if sex == 0 else sex}\n" for index, sex in enumerate(X_SEXES)
            ),
            encoding="utf-8",
        )
    for option in options:
        for threads in (1, 4):
            arguments = ["--bfile", prefix, option]
            if sex_file is not None:
                arguments.extend(("--update-sex", sex_file))
            output = run_legacy(
                gcta,
                Path(work) / option[2:] / f"threads_{threads}",
                "grm",
                arguments,
                threads,
                "Calculating the dominance genetic relationship matrix"
                if kind == "d" else "genetic relationship matrix (GRM) for the X chromosome",
            )
            if option.endswith("gz"):
                assert_gz_grm(output, stem, expected, counts)
            else:
                assert_floats_close(read_floats(path_with(output, f"{stem}.grm.bin")), expected, output)
                assert_floats_close(read_floats(path_with(output, f"{stem}.grm.N.bin")), counts, output)
            outputs[(option, threads)] = output
        if option.endswith("gz"):
            assert_gzip_tables_close(
                path_with(outputs[(option, 1)], f"{stem}.grm.gz"),
                path_with(outputs[(option, 4)], f"{stem}.grm.gz"),
            )
        else:
            for suffix in (f"{stem}.grm.bin", f"{stem}.grm.N.bin", f"{stem}.grm.id"):
                first = path_with(outputs[(option, 1)], suffix)
                second = path_with(outputs[(option, 4)], suffix)
                if suffix.endswith(".id") or first.read_bytes() == second.read_bytes():
                    assert first.read_bytes() == second.read_bytes()
                else:
                    assert_floats_close(read_floats(first), read_floats(second), first)
    if kind == "d":
        for suffix in (".d.grm.bin", ".d.grm.N.bin", ".d.grm.id"):
            assert path_with(outputs[("--make-grm-d-v1", 1)], suffix).read_bytes() == path_with(
                outputs[("--make-grm-d-bin", 1)], suffix
            ).read_bytes()


def scenario_inbred_grm(gcta, work):
    prefix, markers = make_basic_fixture(work)
    standard, counts = expected_legacy_grm(markers, 0)
    expected = [0.5 * value for value in standard]
    outputs = {}
    for option in ("--make-grm-inbred", "--make-grm-inbred-bin", "--make-grm-inbred-gz"):
        for threads in (1, 4):
            output = run_legacy(
                gcta,
                Path(work) / option[2:] / f"threads_{threads}",
                "grm",
                ["--bfile", prefix, option],
                threads,
                "Calculating the genetic relationship matrix",
            )
            if option.endswith("gz"):
                assert_gz_grm(output, "", expected, counts)
            else:
                assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected, output)
                assert_floats_close(read_floats(output.with_suffix(".grm.N.bin")), counts, output)
            outputs[(option, threads)] = output
        if option.endswith("gz"):
            assert_gzip_tables_close(
                path_with(outputs[(option, 1)], ".grm.gz"), path_with(outputs[(option, 4)], ".grm.gz")
            )
        else:
            assert_grm_equal(outputs[(option, 1)], outputs[(option, 4)])
    assert_grm_equal(outputs[("--make-grm-inbred", 1)], outputs[("--make-grm-inbred-bin", 1)])


def expected_recode(markers, mode):
    mus = []
    scales = []
    for marker in markers:
        observed = [value for value in marker if value is not None]
        mu = sum(observed) / len(observed)
        mus.append(mu)
        variance = mu * (1.0 - 0.5 * mu)
        scales.append(0.0 if abs(variance) < 1e-50 else 1.0 / math.sqrt(variance))
    rows = []
    for sample in range(len(markers[0])):
        row = []
        for marker, mu, scale in zip(markers, mus, scales):
            value = marker[sample]
            if mode == "std":
                row.append(0.0 if value is None else (value - mu) * scale)
            elif value is None:
                row.append(mu if mode == "nomiss" else "NA")
            else:
                row.append(float(value))
        rows.append(row)
    return mus, rows


def scenario_recode(gcta, work):
    prefix, markers = make_basic_fixture(work)
    for option, mode in (("--recode", "raw"), ("--recode-nomiss", "nomiss"), ("--recode-std", "std")):
        mus, expected = expected_recode(markers, mode)
        outputs = []
        for threads in (1, 4):
            output = run_legacy(
                gcta,
                Path(work) / option[2:] / f"threads_{threads}",
                "recode",
                ["--bfile", prefix, option],
                threads,
                "Saving the recoded genotype matrix",
            )
            result = path_with(output, ".xmat.gz")
            rows = gzip_table(result)
            assert rows[0] == ["FID", "IID"] + [f"rs{index}" for index in range(len(markers))]
            assert rows[1] == ["Reference", "Allele"] + ["A"] * len(markers)
            for sample, (row, reference) in enumerate(zip(rows[2:], expected)):
                assert row[:2] == ["F1", f"I{sample}"]
                for marker_index, value in enumerate(reference):
                    observed = row[marker_index + 2]
                    if value == "NA":
                        assert observed == "NA"
                    else:
                        assert_close(float(observed), value, f"recode {mode} {sample}/{marker_index}")
            outputs.append(result)
        assert_gzip_tables_close(outputs[0], outputs[1])


def correlation(first, second):
    pairs = [(a, b) for a, b in zip(first, second) if a is not None and b is not None]
    first_mean = sum(a for a, _ in pairs) / len(pairs)
    second_mean = sum(b for _, b in pairs) / len(pairs)
    covariance = sum((a - first_mean) * (b - second_mean) for a, b in pairs)
    first_ss = sum((a - first_mean) ** 2 for a, _ in pairs)
    second_ss = sum((b - second_mean) ** 2 for _, b in pairs)
    if first_ss == 0.0 or second_ss == 0.0:
        return 0.0
    return covariance / math.sqrt(first_ss * second_ss)


def scenario_ld(gcta, work):
    markers = synthetic_markers(40, 20)
    prefix = Path(work) / "input" / "ld"
    write_plink(prefix, markers)
    target = Path(work) / "input" / "target.txt"
    target.write_text("rs5\n", encoding="utf-8")
    outputs = []
    expected_snps = ["rs1", "rs4", "rs7", "rs10", "rs13", "rs16"]
    expected_r = correlation(markers[5], markers[1])
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "ld",
            ["--bfile", prefix, "--ld", target, "--ld-wind", 1000, "--ld-sig", 1.0],
            threads,
            "Estimating LD structure",
        )
        summary = read_table(path_with(output, ".rsq.ld"))
        assert len(summary) == 2
        assert summary[1][0] == "rs5" and summary[1][6] == "6"
        assert_close(float(summary[1][7]), expected_r * expected_r, "mean LD rsq")
        assert_close(float(summary[1][8]), expected_r * expected_r, "median LD rsq")
        assert_close(float(summary[1][9]), expected_r * expected_r, "maximum LD rsq")
        assert read_table(path_with(output, ".snp.ld"))[0] == expected_snps
        for value in read_table(path_with(output, ".r.ld"))[0]:
            assert_close(float(value), expected_r, "signed LD correlation")
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".rsq.ld", ".r.ld", ".snp.ld"))


def scenario_ld_multi(gcta, work):
    prefix = Path(work) / "input" / "ld"
    write_plink(prefix, LD_MARKERS)
    set_groups = ([0, 1], [2, 3])
    filenames = []
    for set_index, group in enumerate(set_groups):
        filename = Path(work) / "input" / f"set{set_index}.txt"
        filename.write_text("".join(f"rs{index}\n" for index in group), encoding="utf-8")
        filenames.append(filename)
    file_list = Path(work) / "input" / "sets.txt"
    file_list.write_text("".join(f"{filename}\n" for filename in filenames), encoding="utf-8")
    rsq = [[correlation(first, second) ** 2 for second in LD_MARKERS] for first in LD_MARKERS]
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "multi",
            ["--bfile", prefix, "--ld-score-multi", file_list, "--ld-wind", 1000],
            threads,
            "Calculating multi-component LD score",
        )
        rows = [row for row in read_table(path_with(output, ".multi_score.ld")) if row]
        assert len(rows) == len(LD_MARKERS) + 1
        for marker_index, row in enumerate(rows[1:]):
            assert row[:3] == [f"rs{marker_index}", "1", str(1000 + marker_index)]
            observed = LD_MARKERS[marker_index]
            frequency = sum(observed) / (2.0 * len(observed))
            assert_close(float(row[3]), min(frequency, 1.0 - frequency), "multi LD MAF")
            for set_index, group in enumerate(set_groups):
                values = [rsq[marker_index][other] for other in group if other != marker_index]
                base = 4 + 3 * set_index
                assert_close(float(row[base]), sum(values) / len(values), "multi LD mean")
                assert_close(float(row[base + 1]), float(len(values)), "multi LD count")
                assert_close(float(row[base + 2]), max(values), "multi LD maximum")
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".multi_score.ld",))


def scenario_ld_region(gcta, work):
    # Multiple chromosomes avoid the legacy single-block regional segmentation dead path.
    markers, chromosomes = family_markers(sample_count=32, marker_count=99)
    prefix = Path(work) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "region",
            ["--bfile", prefix, "--ld-score-region", 10, "--ld-wind", 1000],
            threads,
            "Calculating regional mean LD score",
        )
        rows = read_table(output.with_suffix(".score.ld"))
        assert rows[0] == [
            "SNP", "chr", "bp", "freq", "mean_rsq", "snp_num", "max_rsq", "ldscore_SNP", "ldscore_region"
        ]
        assert len(rows) == len(markers) + 1
        for index, row in enumerate(rows[1:]):
            assert row[:3] == [f"rs{index}", str(chromosomes[index]), str(1000 + index)]
            values = [float(value) for value in row[4:]]
            assert all(math.isfinite(value) for value in values)
            assert_close(values[3], 1.0 + values[0] * values[1], f"regional LD identity {index}")
            assert values[4] > 0.0
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".score.ld",))


def scenario_ld_max(gcta, work):
    markers = synthetic_markers(40, 9)
    prefix = Path(work) / "input" / "ld"
    write_plink(prefix, markers)
    rsq = [[correlation(first, second) ** 2 for second in markers] for first in markers]
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "max",
            ["--bfile", prefix, "--ld-max-rsq", "--ld-wind", 1000],
            threads,
            "Calculating maximum LD rsq",
        )
        rows = [row for row in read_table(path_with(output, ".max_rsq.ld")) if row]
        assert len(rows) == len(markers) + 1
        for index, row in enumerate(rows[1:]):
            assert row[0] == f"rs{index}"
            frequency = sum(markers[index]) / (2.0 * len(markers[index]))
            assert_close(float(row[1]), frequency, "max LD frequency")
            values = [rsq[index][other] for other in range(len(markers)) if other != index]
            if max(values) == 0.0:
                assert float(row[2]) == 0.0
                assert abs(float(row[4])) <= 1e-6
                assert math.isfinite(float(row[5]))
            else:
                assert_close(float(row[2]), max(values), f"max pairwise rsq {index}")
                assert row[3].startswith("rs")
                assert 0.0 <= float(row[4]) <= 1.0
                assert math.isfinite(float(row[5]))
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".max_rsq.ld",))


def scenario_blup(gcta, work):
    prefix, markers = make_polymorphic_fixture(work)
    sample_count = len(markers[0])
    genetic_values = [(index - 3.5) / 3.0 + 0.2 * math.sin(index + 1) for index in range(sample_count)]
    blup_file = Path(work) / "input" / "individual.blup"
    blup_file.write_text(
        "".join(f"F1 I{index} {value:.17g} 0.1\n" for index, value in enumerate(genetic_values)),
        encoding="utf-8",
    )
    expected = []
    for marker in markers:
        observed = [value for value in marker if value is not None]
        mu = sum(observed) / len(observed)
        variance = mu * (1.0 - 0.5 * mu)
        total = sum(
            (value - mu) * genetic_values[sample]
            for sample, value in enumerate(marker) if value is not None
        )
        count = len(observed)
        expected.append(total / variance / count * sample_count / len(markers))
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "blup",
            ["--bfile", prefix, "--blup-snp", blup_file],
            threads,
            "Calculating the BLUP solutions to SNP effects",
        )
        rows = read_table(path_with(output, ".snp.blp"))
        assert len(rows) == len(markers)
        for index, (row, reference) in enumerate(zip(rows, expected)):
            assert row[:2] == [f"rs{index}", "A"]
            assert_close(float(row[2]), reference, f"BLUP SNP {index}")
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".snp.blp",))


def normal_pvalue(z_value):
    return math.erfc(abs(z_value) / math.sqrt(2.0))


def make_summary_fixture(work):
    rng = random.Random(8128)
    sample_count = 48
    marker_count = 12
    markers = []
    for _ in range(marker_count):
        while True:
            marker = [rng.randrange(3) for _ in range(sample_count)]
            if len(set(marker)) == 3:
                break
        markers.append(marker)
    phenotypes = [
        1.3 * (markers[1][sample] - 1.0) - 0.9 * (markers[4][sample] - 1.0) + rng.gauss(0.0, 0.3)
        for sample in range(sample_count)
    ]
    statistics = []
    for marker in markers:
        marker_mean = sum(marker) / sample_count
        phenotype_mean = sum(phenotypes) / sample_count
        variance = sum((value - marker_mean) ** 2 for value in marker)
        beta = sum(
            (marker[sample] - marker_mean) * (phenotypes[sample] - phenotype_mean)
            for sample in range(sample_count)
        ) / variance
        residuals = [
            phenotypes[sample] - phenotype_mean - beta * (marker[sample] - marker_mean)
            for sample in range(sample_count)
        ]
        se = math.sqrt(sum(value * value for value in residuals) / (sample_count - 2) / variance)
        z_value = beta / se
        statistics.append({
            "beta": beta,
            "se": se,
            "p": normal_pvalue(z_value),
            "freq": marker_mean / 2.0,
            "chisq": z_value * z_value,
        })
    meta = Path(work) / "input.ma"
    meta.parent.mkdir(parents=True, exist_ok=True)
    meta.write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G {stat['freq']:.17g} {stat['beta']:.17g} {stat['se']:.17g} "
            f"{stat['p']:.17g} {sample_count}\n"
            for index, stat in enumerate(statistics)
        ),
        encoding="utf-8",
    )
    prefix = Path(work) / "input" / "summary"
    write_plink(prefix, markers)
    return prefix, markers, meta, statistics


def assert_jma(output, statistics, expected_ids):
    rows = read_table(path_with(output, ".jma.cojo"))
    assert len(rows) == len(expected_ids) + 1
    assert [row[1] for row in rows[1:]] == expected_ids
    for row in rows[1:]:
        index = int(row[1][2:])
        assert_close(float(row[4]), statistics[index]["freq"], "JMA frequency")
        assert_close(float(row[5]), statistics[index]["beta"], "JMA marginal beta")
        assert_close(float(row[6]), statistics[index]["se"], "JMA marginal SE")
        assert_close(float(row[7]), statistics[index]["p"], "JMA marginal p")
        assert all(math.isfinite(float(row[column])) for column in (10, 11, 12, 13))
    ld_rows = read_table(path_with(output, ".ldr.cojo"))
    assert len(ld_rows) == len(expected_ids) + 1
    for index in range(len(expected_ids)):
        assert_close(float(ld_rows[index + 1][index + 1]), 1.0, "JMA LD diagonal")


def scenario_cojo_select(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    variants = (
        ("slct", ["--cojo-slct", "--cojo-p", 0.05], ["rs1", "rs4"], "Performing stepwise model selection"),
        ("stepwise", ["--cojo-stepwise", "--cojo-p", 0.05], ["rs1", "rs4"], "Performing stepwise model selection"),
        ("forward", ["--cojo-forward", "--cojo-p", 0.05], ["rs1", "rs4"], "Performing forward model selection"),
        ("backward", ["--cojo-backward", "--cojo-p", 0.05], ["rs1", "rs3", "rs4"], "Performing backward selection"),
        ("top", ["--cojo-top-SNPs", 2], ["rs1", "rs4"], "threshold p-value has been set to 0.5"),
    )
    for name, options, expected_ids, terminal_log in variants:
        outputs = []
        for threads in (1, 4):
            output = run_legacy(
                gcta,
                Path(work) / name / f"threads_{threads}",
                "cojo",
                ["--bfile", prefix, "--cojo-file", meta, "--cojo-collinear", 0.9, "--cojo-wind", 100, "--restrict-output-pC", 1, "--cojo-actual-geno", "--cojo-gc", 1] + options,
                threads,
                terminal_log,
            )
            assert_jma(output, statistics, expected_ids)
            outputs.append(output)
        # The legacy conditional output has uninitialized padding fields after selection.
        compare_suffixes(outputs[0], outputs[1], (".jma.cojo", ".ldr.cojo"))


def scenario_cojo_cond(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    given = Path(work) / "given.txt"
    given.write_text("rs1\n", encoding="utf-8")
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "conditional",
            ["--bfile", prefix, "--cojo-file", meta, "--cojo-cond", given],
            threads,
            "Performing single-SNP association analysis conditional",
        )
        given_rows = read_table(path_with(output, ".given.cojo"))
        assert len(given_rows) == 2 and given_rows[1][1] == "rs1"
        for column, key in ((5, "beta"), (6, "se"), (7, "p")):
            assert_close(float(given_rows[1][column]), statistics[1][key], "given COJO statistic")
        rows = read_table(path_with(output, ".cma.cojo"))
        assert len(rows) == len(statistics)
        assert [row[1] for row in rows[1:]] == [f"rs{index}" for index in range(len(statistics)) if index != 1]
        for row in rows[1:]:
            assert all(math.isfinite(float(row[column])) for column in (10, 11, 12))
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".given.cojo", ".cma.cojo"))


def scenario_cojo_sblup(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "sblup",
            ["--bfile", prefix, "--cojo-file", meta, "--cojo-sblup", 100],
            threads,
            "Estimating the joint effects of all SNPs",
        )
        rows = read_table(path_with(output, ".sblup.cojo"))
        assert len(rows) == len(statistics)
        ratios = []
        for index, row in enumerate(rows):
            assert row[:2] == [f"rs{index}", "A"]
            assert_close(float(row[2]), statistics[index]["beta"], "SBLUP marginal beta")
            assert math.isfinite(float(row[3]))
            ratios.append(float(row[3]))
        assert all(math.isfinite(value) for value in ratios)
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".sblup.cojo",))


def write_snp_sets(work):
    sets = Path(work) / "sets.txt"
    sets.write_text("SET1 rs1 rs2 rs3 END\nSET2 rs4 rs5 END\n", encoding="utf-8")
    return sets


def scenario_fastbat_set(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    sets = write_snp_sets(work)
    groups = (("SET1", (1, 2, 3)), ("SET2", (4, 5)))
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "fastbat",
            ["--bfile", prefix, "--fastBAT", meta, "--fastBAT-set-list", sets, "--fastBAT-write-snpset", "--fastBAT-ld-cutoff", 0.9, "--fastBAT-wind", 50],
            threads,
            "Running fastBAT analysis",
        )
        rows = read_table(path_with(output, ".fastbat"))
        assert len(rows) == 3
        for row, (name, indices) in zip(rows[1:], groups):
            assert row[0] == name and row[1] == str(len(indices))
            assert_close(float(row[2]), sum(statistics[index]["chisq"] for index in indices), "fastBAT chisq")
            best = min(indices, key=lambda index: statistics[index]["p"])
            assert_close(float(row[4]), statistics[best]["p"], "fastBAT top p")
            assert row[5] == f"rs{best}"
            assert 0.0 < float(row[3]) < 1.0
        assert path_with(output, ".snpset").read_text(encoding="utf-8") == (
            "SET1\nrs1\nrs2\nrs3\nEND\n\nSET2\nrs4\nrs5\nEND\n\n"
        )
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".fastbat", ".snpset"))


def scenario_fastbat_seg(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "segment",
            ["--bfile", prefix, "--fastBAT", meta, "--fastBAT-seg", 10, "--fastBAT-write-snpset"],
            threads,
            "Running fastBAT analysis at genomic segments",
        )
        rows = read_table(path_with(output, ".seg.fastbat"))
        assert len(rows) == 2
        assert rows[1][:4] == ["1", "1000", "1011", str(len(statistics))]
        assert_close(float(rows[1][4]), sum(stat["chisq"] for stat in statistics), "segment chisq")
        best = min(range(len(statistics)), key=lambda index: statistics[index]["p"])
        assert_close(float(rows[1][6]), statistics[best]["p"], "segment top p")
        assert rows[1][7] == f"rs{best}"
        assert 0.0 < float(rows[1][5]) < 1.0
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".seg.fastbat", ".seg.snpset"))


def scenario_mbat_set(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    sets = write_snp_sets(work)
    groups = (("SET1", (1, 2, 3)), ("SET2", (4, 5)))
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "mbat",
            ["--bfile", prefix, "--mBAT-combo", meta, "--mBAT-set-list", sets, "--mBAT-write-snpset", "--mBAT-svd-gamma", 0.9, "--mBAT-wind", 50, "--mBAT-print-all-p"],
            threads,
            "Running mBAT-combo analysis",
        )
        rows = read_table(path_with(output, ".set.assoc.mbat"))
        assert len(rows) == 3
        for row, (name, indices) in zip(rows[1:], groups):
            assert row[0] == name and row[1] == str(len(indices))
            best = min(indices, key=lambda index: statistics[index]["p"])
            assert row[2] == f"rs{best}"
            assert_close(float(row[3]), statistics[best]["p"], "mBAT top p")
            assert 1 <= int(row[4]) <= len(indices)
            assert float(row[5]) > 0.0 and 0.0 < float(row[6]) < 1.0
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".set.assoc.mbat", ".set.snpset.mbat"))


def write_eigen_files(work, sample_count):
    prefix = Path(work) / "pcs"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".eigenval").write_text("2\n1\n", encoding="utf-8")
    prefix.with_suffix(".eigenvec").write_text(
        "".join(
            f"F1 I{index} {math.cos(index + 1):.17g} {math.sin(index + 1):.17g}\n"
            for index in range(sample_count)
        ),
        encoding="utf-8",
    )
    return prefix


def scenario_gwas_adj_pc(gcta, work):
    prefix, _, meta, statistics = make_summary_fixture(work)
    pcs = write_eigen_files(work, 48)
    covariate = Path(work) / "covariate.ma"
    covariate_stats = []
    for index, stat in enumerate(statistics):
        beta = 0.35 * stat["beta"] + 0.02 * (index % 3 - 1)
        se = stat["se"]
        covariate_stats.append((beta, se, normal_pvalue(beta / se)))
    covariate.write_text(
        "SNP A1 A2 freq b se p N\n"
        + "".join(
            f"rs{index} A G {statistics[index]['freq']:.17g} {beta:.17g} {se:.17g} {p:.17g} 48\n"
            for index, (beta, se, p) in enumerate(covariate_stats)
        ),
        encoding="utf-8",
    )
    file_list = Path(work) / "pcadjust.txt"
    file_list.write_text(f"{meta}\n{covariate}\n", encoding="utf-8")
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "pcadj",
            ["--bfile", prefix, "--gwas-adj-pc", file_list, "--gwas-adj-pc-wind", 10000, "--pc-loading", pcs],
            threads,
            "PC1, bxy =",
        )
        rows = read_table(path_with(output, ".pcadj.cma"))
        assert len(rows) == len(statistics) + 1
        changed = False
        for index, row in enumerate(rows[1:]):
            assert row[0] == f"rs{index}"
            assert_close(float(row[3]), statistics[index]["freq"], "PC-adjust frequency")
            assert_close(float(row[4]), statistics[index]["beta"], "PC-adjust original beta")
            assert_close(float(row[5]), statistics[index]["se"], "PC-adjust original SE")
            assert_close(float(row[6]), statistics[index]["p"], "PC-adjust original p")
            assert math.isfinite(float(row[7]))
            assert math.isfinite(float(row[8])) and math.isfinite(float(row[9]))
            changed |= not math.isclose(float(row[8]), statistics[index]["beta"], rel_tol=1e-4, abs_tol=1e-6)
        assert changed
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".pcadj.cma",))


def scenario_pc_loading(gcta, work):
    prefix, markers = make_polymorphic_fixture(work)
    pcs = write_eigen_files(work, len(markers[0]))
    eigenvectors = [(math.cos(index + 1), math.sin(index + 1)) for index in range(len(markers[0]))]
    eigenvalues = (2.0, 1.0)
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "loading",
            ["--bfile", prefix, "--pc-loading", pcs],
            threads,
            "Calculating PC loadings of SNPs",
        )
        rows = read_table(path_with(output, ".pcl"))
        assert rows[0] == ["SNP", "A1", "A2", "mu", "pc1_loading", "pc2_loading"]
        for marker_index, (marker, row) in enumerate(zip(markers, rows[1:])):
            observed = [value for value in marker if value is not None]
            mu = sum(observed) / len(observed)
            scale = 1.0 / math.sqrt(mu * (1.0 - 0.5 * mu))
            standardized = [0.0 if value is None else (value - mu) * scale for value in marker]
            assert row[:3] == [f"rs{marker_index}", "A", "G"]
            assert_close(float(row[3]), mu, "PC loading mu")
            for pc in range(2):
                expected = sum(
                    eigenvectors[sample][pc] * standardized[sample]
                    for sample in range(len(marker))
                ) / (eigenvalues[pc] * len(markers))
                assert_close(float(row[4 + pc]), expected, f"PC loading {marker_index}/{pc}")
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".pcl",))


def scenario_project_loading(gcta, work):
    prefix, markers = make_polymorphic_fixture(work)
    mus = []
    for marker in markers:
        observed = [value for value in marker if value is not None]
        mus.append(sum(observed) / len(observed))
    loadings = [(0.03 * (index + 1), -0.02 + 0.005 * index) for index in range(len(markers))]
    weights = Path(work) / "weights.pcl"
    weights.write_text(
        "SNP\tA1\tA2\tmu\tpc1_loading\tpc2_loading\n"
        + "".join(
            f"rs{index}\tA\tG\t{mus[index]:.17g}\t{loadings[index][0]:.17g}\t{loadings[index][1]:.17g}\n"
            for index in range(len(markers))
        )
        + "rsMissing\tA\tG\t1\t0.1\t0.1\n",
        encoding="utf-8",
    )
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "project",
            ["--bfile", prefix, "--project-loading", weights.with_suffix(""), 2],
            threads,
            "Standardizing genotypes and projecting PCs",
        )
        assert path_with(output, ".proj.missnp").read_text(encoding="utf-8") == "rsMissing\n"
        rows = read_table(path_with(output, ".proj.eigenvec"))
        assert len(rows) == len(markers[0])
        for sample, row in enumerate(rows):
            assert row[:2] == ["F1", f"I{sample}"]
            for pc in range(2):
                expected = 0.0
                for marker, mu, loading in zip(markers, mus, loadings):
                    value = marker[sample]
                    if value is not None:
                        expected += (value - mu) / math.sqrt(mu * (1.0 - 0.5 * mu)) * loading[pc]
                assert_close(float(row[2 + pc]), expected, f"projection {sample}/{pc}")
        outputs.append(output)
    compare_suffixes(outputs[0], outputs[1], (".proj.eigenvec", ".proj.missnp"))


def make_family_grm(gcta, work, sample_count=64, marker_count=201):
    markers, chromosomes = family_markers(sample_count=sample_count, marker_count=marker_count)
    prefix = Path(work) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, _ = write_family_phenotypes(prefix, markers)
    grm = run_gcta(gcta, work, "grm", ["--bfile", prefix, "--make-grm-v1"])
    return grm, pheno


def scenario_hereg_bivar(gcta, work):
    grm, pheno = make_family_grm(gcta, work, sample_count=32, marker_count=99)
    two_traits = Path(work) / "two.phen"
    rows = []
    for index, line in enumerate(pheno.read_text(encoding="utf-8").splitlines()):
        fields = line.split()
        first = float(fields[2])
        rows.append(f"{fields[0]} {fields[1]} {first:.17g} {0.8 * first + 0.05 * (index % 7 - 3):.17g}\n")
    two_traits.write_text("".join(rows), encoding="utf-8")
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "he",
            ["--grm", grm, "--pheno", two_traits, "--HEreg-bivar", 1, 2],
            threads,
            "Performing Haseman-Elston regression",
        )
        result = path_with(output, ".HEreg")
        text = result.read_text(encoding="utf-8")
        assert "HE-CP" in text and "V(G)/Vp_tr1" in text and "C(G)/Vp_tr12" in text
        values = {}
        for row in read_table(result):
            if len(row) >= 2:
                try:
                    values[row[0]] = float(row[1])
                except ValueError:
                    pass
        assert values["N_tr1"] == values["N_tr2"] == 32
        for name in ("V(G)/Vp_tr1", "V(G)/Vp_tr2", "C(G)/Vp_tr12"):
            assert math.isfinite(values[name])
        outputs.append(result)
    assert_tables_close(outputs[0], outputs[1])


def assert_hsq(path, sample_count):
    values = {}
    for row in read_table(path):
        if len(row) < 2:
            continue
        try:
            values[row[0]] = float(row[1])
        except ValueError:
            pass
    assert values["n"] == sample_count
    for name in ("V(G)", "V(e)", "Vp", "V(G)/Vp", "logL"):
        assert math.isfinite(values[name]), (path, name, values)
    assert_close(values["V(G)"] + values["V(e)"], values["Vp"], "REML variance sum")
    assert_close(values["V(G)"] / values["Vp"], values["V(G)/Vp"], "REML heritability")


def scenario_reml_variants(gcta, work):
    grm, pheno = make_family_grm(gcta, work)
    for option in ("--reml-no-constrain", "--reml-wfam"):
        outputs = []
        for threads in (1, 4):
            output = run_legacy(
                gcta,
                Path(work) / option[2:] / f"threads_{threads}",
                "reml",
                [
                    "--grm", grm, "--pheno", pheno, option, "--reml-priors", 0.5,
                    "--reml-no-lrt", "--reml-maxit", 100,
                ],
                threads,
                "Log-likelihood ratio converged.",
            )
            assert_hsq(output.with_suffix(".hsq"), 64)
            outputs.append(output)
        compare_suffixes(outputs[0], outputs[1], (".hsq",))


def scenario_grm_align(gcta, work):
    prefix, _ = make_basic_fixture(work)
    full = run_gcta(gcta, work, "full", ["--bfile", prefix, "--make-grm-v1"])
    keep = Path(work) / "keep.txt"
    keep.write_text("".join(f"F1 I{index}\n" for index in range(2, 8)), encoding="utf-8")
    subset = run_gcta(gcta, work, "subset", ["--bfile", prefix, "--keep", keep, "--make-grm-v1"])
    grm_list = Path(work) / "grms.txt"
    grm_list.write_text(f"{full}\n{subset}\n", encoding="utf-8")
    for threads in (1, 4):
        case = Path(work) / f"threads_{threads}"
        case.mkdir(parents=True, exist_ok=True)
        run_legacy(
            gcta,
            case,
            "align",
            ["--mgrm", grm_list, "--grm-align"],
            threads,
            "2 GRMs have been aligned",
        )
        assert path_with(full, ".aligned.grm.id").read_bytes() == path_with(subset, ".aligned.grm.id").read_bytes()
        assert all(math.isfinite(value) for value in read_floats(path_with(full, ".aligned.grm.bin")))
        assert all(math.isfinite(value) for value in read_floats(path_with(subset, ".aligned.grm.bin")))


def scenario_make_bk(gcta, work):
    prefix, _ = make_basic_fixture(work)
    grm = run_gcta(gcta, work, "source", ["--bfile", prefix, "--make-grm-v1"])
    original = read_floats(path_with(grm, ".grm.bin"))
    expected = []
    sample_count = 8
    index = 0
    for row in range(sample_count):
        for column in range(row + 1):
            value = original[index]
            expected.append(0.0 if column < row and value < 0.05 else value)
            index += 1
    outputs = []
    for threads in (1, 4):
        output = run_legacy(
            gcta,
            Path(work) / f"threads_{threads}",
            "bk",
            ["--grm", grm, "--make-bK", 0.05, "--save-ram"],
            threads,
            "Pruning the GRM to a sparse matrix with a cutoff of 0.050000",
        )
        assert_floats_close(read_floats(output.with_suffix(".grm.bin")), expected, output)
        assert path_with(output, ".grm.N.bin").read_bytes() == path_with(grm, ".grm.N.bin").read_bytes()
        assert path_with(output, ".grm.id").read_bytes() == path_with(grm, ".grm.id").read_bytes()
        outputs.append(output)
    assert_grm_equal(outputs[0], outputs[1])


def main():
    parser = argparse.ArgumentParser(description="Behavioral tests for legacy terminal commands")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    args = parser.parse_args()
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)

    scenario = args.scenario
    if scenario == "raw-files":
        scenario_raw_files(args.gcta, work)
    elif scenario == "hapmap-genet-dst":
        scenario_hapmap(args.gcta, work)
    elif scenario == "paa":
        scenario_paa(args.gcta, work)
    elif scenario == "ibc":
        scenario_ibc(args.gcta, work)
    elif scenario == "make-grm-formats":
        scenario_make_grm_formats(args.gcta, work)
    elif scenario == "make-grm-d-formats":
        scenario_special_grm(args.gcta, work, "d")
    elif scenario == "make-grm-xchr-formats":
        scenario_special_grm(args.gcta, work, "xchr")
    elif scenario == "make-grm-inbred-formats":
        scenario_inbred_grm(args.gcta, work)
    elif scenario == "recode":
        scenario_recode(args.gcta, work)
    elif scenario == "ld":
        scenario_ld(args.gcta, work)
    elif scenario == "ld-score-multi":
        scenario_ld_multi(args.gcta, work)
    elif scenario == "ld-score-region":
        scenario_ld_region(args.gcta, work)
    elif scenario == "ld-max-rsq":
        scenario_ld_max(args.gcta, work)
    elif scenario == "blup-snp":
        scenario_blup(args.gcta, work)
    elif scenario == "cojo-select":
        scenario_cojo_select(args.gcta, work)
    elif scenario == "cojo-cond":
        scenario_cojo_cond(args.gcta, work)
    elif scenario == "cojo-sblup":
        scenario_cojo_sblup(args.gcta, work)
    elif scenario == "fastbat-set":
        scenario_fastbat_set(args.gcta, work)
    elif scenario == "fastbat-seg":
        scenario_fastbat_seg(args.gcta, work)
    elif scenario == "mbat-set":
        scenario_mbat_set(args.gcta, work)
    elif scenario == "gwas-adj-pc":
        scenario_gwas_adj_pc(args.gcta, work)
    elif scenario == "pc-loading":
        scenario_pc_loading(args.gcta, work)
    elif scenario == "project-loading":
        scenario_project_loading(args.gcta, work)
    elif scenario == "hereg-bivar":
        scenario_hereg_bivar(args.gcta, work)
    elif scenario == "reml-variants":
        scenario_reml_variants(args.gcta, work)
    elif scenario == "grm-align":
        scenario_grm_align(args.gcta, work)
    else:
        scenario_make_bk(args.gcta, work)
    print(f"Passed legacy terminal regression scenario: {scenario}")


if __name__ == "__main__":
    main()
