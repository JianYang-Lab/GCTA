#!/usr/bin/env python3

"""Behavioral coverage for v2 options not covered by the focused suites.

Most numeric checks use independently calculated values.  A few options are
currently allow-listed but have no active v2 implementation (``--bed``,
``--dosage``, ``--noblas``, ``--num-vec``, and ``--trace-exact``).  The
``allowlisted_noops`` scenario explicitly characterizes that state: it checks
v2 dispatch and either the exact missing-loader diagnostic or equivalence to
the same calculation without the no-op.  It does not claim those options have
an implemented numerical effect.  Likewise, ``--r2`` currently selects the
r2 method flag but LD.cpp does not consume it; the LD scenario records the
observed equivalence to r rather than silently treating r2 as implemented.
"""

import argparse
import math
from pathlib import Path
import struct

from bgen import write_bgen
from gcta_runner import (
    assert_floats_close,
    assert_grm_equal,
    assert_tables_close,
    read_floats,
    read_table,
    run_gcta,
)
from pgen import write_pgen
from plink import expected_grm, synthetic_markers, write_plink
from synthetic import family_markers, write_family_phenotypes
from test_acat import expected_acat, write_inputs as write_acat_inputs
from test_errors import expect_error
from test_expression import write_ecojo_fixture
from test_grm_x import MARKERS as X_MARKERS, SEXES as X_SEXES
from test_gsmr import make_fixture as make_gsmr_fixture
from test_ld import expected_values as expected_ld_values, read_zld


SCENARIOS = (
    "components",
    "filters",
    "pruning",
    "partitions",
    "multifile",
    "pgen_inputs",
    "bgen_dosage",
    "grm_merge",
    "ld",
    "fastgwa_lr",
    "fastgwa_exact",
    "fastgwa_grammar",
    "fastgwa_binary",
    "fastgwa_ge",
    "inversion_errors",
    "obsolete_errors",
    "burden",
    "gsmr",
    "gsmr_errors",
    "ecojo",
    "acat_wind",
    "allowlisted_noops",
)


def log_text(output):
    return Path(f"{output}.log").read_text(encoding="utf-8")


def console_text(output):
    return Path(f"{output}.console").read_text(encoding="utf-8")


def assert_v2(output, marker):
    log = log_text(output)
    assert "Options:" in log, f"Missing v2 option listing: {output}.log"
    assert marker in log, f"Missing v2 marker {marker!r}: {output}.log"
    assert "Accepted options:" not in console_text(output), f"Legacy fallback: {output}"


def assert_option_listed(output, *flags):
    log = log_text(output)
    for flag in flags:
        assert f"\n{flag}" in log, f"Option was not dispatched through v2: {flag}"


def expect_v2_error(gcta, work, name, arguments, diagnostic):
    work = Path(work)
    output = work / name
    expect_error(gcta, work, name, [*arguments, "--out", output])
    console = (work / f"{name}.console").read_text(encoding="utf-8")
    assert diagnostic in console, f"Missing diagnostic {diagnostic!r}:\n{console}"
    assert "Accepted options:" not in console, f"Unexpected legacy fallback: {name}"


def write_ids(path, ids):
    Path(path).write_text("".join(f"F1\tI{index}\n" for index in ids), encoding="utf-8")


def write_manual_grm(prefix, ids, matrix, counts=None):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    write_ids(Path(f"{prefix}.grm.id"), ids)
    values = []
    count_values = []
    for row in range(len(ids)):
        for column in range(row + 1):
            values.append(matrix[row][column])
            count_values.append(10.0 if counts is None else counts[row][column])
    Path(f"{prefix}.grm.bin").write_bytes(struct.pack(f"<{len(values)}f", *values))
    Path(f"{prefix}.grm.N.bin").write_bytes(struct.pack(f"<{len(values)}f", *count_values))


def read_ids(path):
    return [int(line.split()[1][1:]) for line in Path(path).read_text(encoding="utf-8").splitlines()]


def run_components(gcta, work, threads):
    work = Path(work)
    markers = synthetic_markers(sample_count=8, marker_count=12)
    prefix = work / "input" / "components"
    write_plink(prefix, markers)
    grm = run_gcta(
        gcta,
        work,
        "grm",
        [
            "--bfile", prefix,
            "--bim", prefix.with_suffix(".bim"),
            "--fam", prefix.with_suffix(".fam"),
            "--make-grm",
        ],
        threads=threads,
    )
    assert_v2(grm, "The GRM computation is completed.")
    assert_option_listed(grm, "--bim", "--fam")
    expected, counts = expected_grm(markers, 0)
    assert_floats_close(read_floats(Path(f"{grm}.grm.bin")), expected, grm)
    assert read_floats(Path(f"{grm}.grm.N.bin")) == counts

    references = work / "references.txt"
    references.write_text("".join(f"rs{index} G\n" for index in range(len(markers))), encoding="utf-8")
    covariates = work / "covar.txt"
    covariates.write_text(
        "".join(f"F1 I{index} {'A' if index % 2 else 'B'}\n" for index in range(8)),
        encoding="utf-8",
    )
    updated = run_gcta(
        gcta,
        work,
        "updated",
        [
            "--bfile", prefix,
            "--update-ref-allele", references,
            "--freq",
        ],
        threads=threads,
    )
    assert_v2(updated, "Computing allele frequencies")
    rows = read_table(updated.with_suffix(".frq"))
    assert len(rows) == len(markers) + 1
    for index, row in enumerate(rows[1:]):
        observed = [value for value in markers[index] if value is not None]
        expected_frequency = 1.0 - sum(observed) / (2.0 * len(observed))
        assert row[3:5] == ["G", "A"]
        assert math.isclose(float(row[5]), expected_frequency, rel_tol=2e-5, abs_tol=1e-6)

    frequencies = work / "frequencies.txt"
    supplied = [0.2 + 0.01 * (index % 7) for index in range(len(markers))]
    frequencies.write_text(
        "".join(f"rs{index} A {frequency}\n" for index, frequency in enumerate(supplied)),
        encoding="utf-8",
    )
    supplied_output = run_gcta(
        gcta,
        work,
        "supplied",
        ["--bfile", prefix, "--update-freq", frequencies, "--freq"],
        threads=threads,
    )
    assert_v2(supplied_output, f"Frequencies of {len(markers)} SNPs are updated.")
    rows = read_table(supplied_output.with_suffix(".frq"))
    for row, expected_frequency in zip(rows[1:], supplied):
        assert math.isclose(float(row[5]), expected_frequency, rel_tol=2e-5, abs_tol=1e-6)
    return grm, updated.with_suffix(".frq"), supplied_output.with_suffix(".frq")


def run_filters(gcta, work, threads):
    work = Path(work)
    markers = synthetic_markers(sample_count=8, marker_count=17)
    markers[2][7] = None
    markers[3][2] = None
    markers[4][1] = markers[4][2] = None
    markers[5][0] = markers[5][1] = markers[5][2] = None
    prefix = work / "input" / "filters"
    write_bgen(prefix, markers)
    output = run_gcta(
        gcta,
        work,
        "filtered",
        [
            "--bgen", prefix.with_suffix(".bgen"),
            "--sample", prefix.with_suffix(".sample"),
            "--geno", 0.75,
            "--info", 0.8,
            "--freq",
        ],
        threads=threads,
    )
    assert_v2(output, "Computing allele frequencies")
    assert_option_listed(output, "--geno", "--info")
    rows = read_table(output.with_suffix(".frq"))
    expected_indices = [
        index for index, marker in enumerate(markers)
        if sum(value is None for value in marker) / len(marker) <= 0.75
    ]
    assert [int(row[1][2:]) for row in rows[1:]] == expected_indices
    for row, index in zip(rows[1:], expected_indices):
        observed = [value for value in markers[index] if value is not None]
        assert math.isclose(
            float(row[5]), sum(observed) / (2.0 * len(observed)), rel_tol=2e-5, abs_tol=1e-6
        )
        assert float(row[7]) >= 0.8
    return output.with_suffix(".frq")


def expected_pruned_ids(ids, matrix, threshold):
    edges = [
        (row, column)
        for row in range(len(ids))
        for column in range(row)
        if matrix[row][column] > threshold
    ]
    degrees = {index: 0 for index in range(len(ids))}
    for row, column in edges:
        degrees[row] += 1
        degrees[column] += 1
    removed = set()
    for row, column in edges:
        if degrees[row] < degrees[column]:
            row, column = column, row
        removed.add(row)
    return [index for index in ids if index not in removed], edges


def run_pruning(gcta, work, threads):
    work = Path(work)
    ids = list(range(5))
    matrix = [[0.0] * 5 for _ in range(5)]
    for index in ids:
        matrix[index][index] = 1.0
    for row, column, value in ((1, 0, 0.8), (3, 2, 0.7), (4, 0, 0.6)):
        matrix[row][column] = matrix[column][row] = value
    source = work / "input" / "relationships"
    write_manual_grm(source, ids, matrix)
    kept, edges = expected_pruned_ids(ids, matrix, 0.5)

    cutoff = run_gcta(
        gcta,
        work,
        "cutoff",
        ["--grm", source, "--make-grm", "--grm-cutoff", 0.5, "--cutoff-detail"],
        threads=threads,
    )
    assert_v2(cutoff, "Pruning the GRM with a cutoff")
    assert read_ids(Path(f"{cutoff}.grm.id")) == kept
    assert read_ids(Path(f"{cutoff}.singleton.txt")) == kept
    family_rows = read_table(Path(f"{cutoff}.family.txt"))
    assert len(family_rows) == len(edges)
    for row, (expected_row, expected_column) in zip(family_rows, edges):
        assert [int(row[1][1:]), int(row[3][1:])] == [expected_row, expected_column]

    singleton = run_gcta(
        gcta,
        work,
        "singleton",
        ["--grm", source, "--grm-singleton", 0.5],
        threads=threads,
    )
    assert_v2(singleton, "Pruned singleton IDs")
    assert read_ids(Path(f"{singleton}.singleton.txt")) == kept
    assert not Path(f"{singleton}.grm.bin").exists()

    bk = run_gcta(
        gcta,
        work,
        "bk",
        ["--grm", source, "--make-bK", 0.5],
        threads=threads,
    )
    assert_v2(bk, "GRM has been saved")
    expected_bk = []
    expected_n = []
    for row in ids:
        for column in range(row + 1):
            expected_bk.append(matrix[row][column] if matrix[row][column] > 0.5 else 0.0)
            expected_n.append(10.0)
    assert_floats_close(read_floats(Path(f"{bk}.grm.bin")), expected_bk, bk)
    assert read_floats(Path(f"{bk}.grm.N.bin")) == expected_n
    return cutoff, bk


def make_x_fixture(work):
    prefix = Path(work) / "input" / "x"
    write_plink(prefix, X_MARKERS, chromosomes=[23] * len(X_MARKERS))
    prefix.with_suffix(".fam").write_text(
        "".join(f"F1 I{index} 0 0 0 -9\n" for index in range(len(X_SEXES))), encoding="utf-8"
    )
    sex = Path(work) / "sex.txt"
    sex.write_text(
        "".join(f"F1 I{index} {value}\n" for index, value in enumerate(X_SEXES)), encoding="utf-8"
    )
    return prefix, sex


def run_x_parts(gcta, work, threads):
    prefix, sex = make_x_fixture(work)
    common = ["--bfile", prefix, "--update-sex", sex, "--dc", 0]
    whole = run_gcta(
        gcta, work, "whole", [*common, "--make-grm-xchr"], threads=threads
    )
    assert_v2(whole, "Note: this function takes X chromosome as non-PAR region.")
    assert f"{len(X_SEXES)} individuals with valid sex information to be included" in log_text(whole)
    parts = []
    for part in (1, 2):
        output = run_gcta(
            gcta,
            work,
            f"part_{part}",
            [*common, "--make-grm-xchr-part", 2, part],
            threads=threads,
        )
        assert "Subset" in Path(f"{output}.part_2_{part}.log").read_text(encoding="utf-8")
        assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
        parts.append(Path(f"{output}.part_2_{part}"))
    joined = Path(work) / "joined"
    for suffix in ("grm.bin", "grm.N.bin", "grm.id"):
        Path(f"{joined}.{suffix}").write_bytes(
            b"".join(Path(f"{part}.{suffix}").read_bytes() for part in parts)
        )
    assert_grm_equal(whole, joined)
    return whole


def run_dominance_parts(gcta, work, threads):
    work = Path(work)
    markers = synthetic_markers(sample_count=8, marker_count=17)
    prefix = work / "input" / "dominance"
    write_plink(prefix, markers)
    whole = run_gcta(gcta, work, "whole", ["--bfile", prefix, "--make-grm-d"], threads=threads)
    parts = []
    for part in (1, 2):
        output = run_gcta(
            gcta,
            work,
            f"part_{part}",
            ["--bfile", prefix, "--make-grm-d-part", 2, part],
            threads=threads,
        )
        assert "Subset" in Path(f"{output}.part_2_{part}.log").read_text(encoding="utf-8")
        assert "Accepted options:" not in Path(f"{output}.console").read_text(encoding="utf-8")
        parts.append(Path(f"{output}.part_2_{part}"))
    joined = Path(work) / "joined"
    for suffix in ("d.grm.bin", "d.grm.N.bin", "d.grm.id"):
        Path(f"{joined}.{suffix}").write_bytes(
            b"".join(Path(f"{part}.{suffix}").read_bytes() for part in parts)
        )
    for suffix in ("d.grm.bin", "d.grm.N.bin"):
        assert_floats_close(
            read_floats(Path(f"{whole}.{suffix}")), read_floats(Path(f"{joined}.{suffix}")), joined
        )
    assert Path(f"{whole}.d.grm.id").read_bytes() == Path(f"{joined}.d.grm.id").read_bytes()
    return Path(f"{whole}.d.grm.bin")


def write_multifile_fixture(work, format_name):
    work = Path(work)
    markers = synthetic_markers(sample_count=8, marker_count=12)
    partitions = (markers[:4], markers[4:])
    prefixes = []
    for index, partition in enumerate(partitions):
        prefix = work / "input" / f"part{index}"
        if format_name == "bgen":
            write_bgen(prefix, partition)
        else:
            write_pgen(prefix, partition)
        prefixes.append(prefix)
    file_list = work / "files.txt"
    paths = [path.with_suffix(".bgen") for path in prefixes] if format_name == "bgen" else prefixes
    file_list.write_text("".join(f"{path}\n" for path in paths), encoding="utf-8")
    return markers, prefixes, file_list


def run_multifile(gcta, work, threads):
    markers, _, file_list = write_multifile_fixture(work, "pgen")
    plink_work = Path(work) / "plink"
    plink_prefix = plink_work / "input" / "source"
    write_plink(plink_prefix, markers)
    # --dosage is allow-listed but has no v2 reader; with PLINK input it is a no-op.
    output = run_gcta(
        gcta,
        plink_work,
        "grm",
        ["--mbfile", file_list, "--dosage", "--make-grm"],
        threads=threads,
    )
    assert_v2(output, "The GRM computation is completed.")
    expected, counts = expected_grm(markers, 0)
    assert_floats_close(read_floats(Path(f"{output}.grm.bin")), expected, output)
    assert read_floats(Path(f"{output}.grm.N.bin")) == counts
    return output


def run_pgen_inputs(gcta, work, threads):
    outputs = []
    for option in ("--bpfile", "--mpfile", "--mbpfile"):
        option_work = Path(work) / option[2:]
        markers, prefixes, file_list = write_multifile_fixture(option_work, "pgen")
        source = prefixes[0] if option == "--bpfile" else file_list
        selected_markers = markers[:4] if option == "--bpfile" else markers
        output = run_gcta(
            gcta, option_work, "grm", [option, source, "--make-grm"], threads=threads
        )
        assert_v2(output, "The GRM computation is completed.")
        expected, counts = expected_grm(selected_markers, 0)
        assert_floats_close(read_floats(Path(f"{output}.grm.bin")), expected, output)
        assert read_floats(Path(f"{output}.grm.N.bin")) == counts
        outputs.append(output)
    return outputs


def run_bgen_dosage(gcta, work, threads):
    markers, prefixes, file_list = write_multifile_fixture(work, "bgen")
    output = run_gcta(
        gcta,
        work,
        "grm",
        [
            "--mbgen", file_list,
            "--sample", prefixes[0].with_suffix(".sample"),
            "--hard-call-thresh", 0.9,
            "--dosage-call",
            "--make-grm",
        ],
        threads=threads,
    )
    assert_v2(output, "The GRM computation is completed.")
    expected, counts = expected_grm(markers, 0)
    assert_floats_close(read_floats(Path(f"{output}.grm.bin")), expected, output)
    assert read_floats(Path(f"{output}.grm.N.bin")) == counts
    return output


def run_grm_merge(gcta, work, threads):
    work = Path(work)
    first_ids = [0, 1, 2, 3]
    second_ids = [1, 2, 3, 4]
    first = [[1.0 if row == column else (row + column + 1) / 20.0 for column in range(4)] for row in range(4)]
    second = [[1.0 if row == column else (row * 4 + column + 2) / 25.0 for column in range(4)] for row in range(4)]
    first_prefix = work / "input" / "first"
    second_prefix = work / "input" / "second"
    write_manual_grm(first_prefix, first_ids, first)
    write_manual_grm(second_prefix, second_ids, second)
    file_list = work / "unify.txt"
    file_list.write_text(f"{first_prefix}\n{second_prefix}\n", encoding="utf-8")
    unified = run_gcta(gcta, work, "unified", ["--mgrm", file_list, "--unify-grm"], threads=threads)
    assert_v2(unified, "Writing unified GRM in binary format")
    common = [1, 2, 3]
    outputs = [work / f"{name}_unified" for name in ("first", "second")]
    for output in outputs:
        assert read_ids(Path(f"{output}.grm.id")) == common
    expected_first = []
    expected_second = []
    first_indices = [first_ids.index(index) for index in common]
    second_indices = [second_ids.index(index) for index in common]
    for values, indices, destination in (
        (first, first_indices, expected_first),
        (second, second_indices, expected_second),
    ):
        for row in indices:
            for column in indices:
                if indices.index(column) <= indices.index(row):
                    destination.append(values[row][column])
    assert_floats_close(read_floats(Path(f"{outputs[0]}.grm.bin")), expected_first, outputs[0])
    assert_floats_close(read_floats(Path(f"{outputs[1]}.grm.bin")), expected_second, outputs[1])

    subtract_ids = [0, 1, 2]
    subtract_matrix = [[1.0 if row == column else 0.1 for column in range(3)] for row in range(3)]
    subtract_prefix = work / "input" / "subtract"
    write_manual_grm(subtract_prefix, subtract_ids, subtract_matrix)
    base_prefix = work / "input" / "base"
    base_matrix = [[1.0 if row == column else 0.3 for column in range(3)] for row in range(3)]
    base_counts = [[20.0 for _ in range(3)] for _ in range(3)]
    write_manual_grm(base_prefix, subtract_ids, base_matrix, base_counts)
    subtract_list = work / "subtract.txt"
    subtract_list.write_text(f"{base_prefix}\n{subtract_prefix}\n", encoding="utf-8")
    result = run_gcta(
        gcta, work, "difference", ["--mgrm", subtract_list, "--subtract-grm"], threads=threads
    )
    assert_v2(result, "Subtracting GRMs")
    expected_values = []
    for row in range(3):
        for column in range(row + 1):
            expected_values.append((base_matrix[row][column] * 20.0 - subtract_matrix[row][column] * 10.0) / 10.0)
    assert_floats_close(read_floats(Path(f"{result}.grm.bin")), expected_values, result)
    assert read_floats(Path(f"{result}.grm.N.bin")) == [10.0] * 6
    return outputs[0], result


def run_ld(gcta, work, threads):
    work = Path(work)
    markers = synthetic_markers(sample_count=8, marker_count=12)
    prefix = work / "input" / "ld"
    write_plink(prefix, markers)
    outputs = {}
    expected_r = expected_ld_values(markers, 1000)
    for method in ("r", "r2"):
        output = run_gcta(
            gcta,
            work,
            method,
            ["--bfile", prefix, "--ld-matrix", "--ld-wind", 1000, f"--{method}"],
            threads=threads,
        )
        assert_v2(output, "Generating LD matrix")
        actual = read_zld(output.with_suffix(".zld"))
        for actual_row, expected_row in zip(actual, expected_r):
            # Characterization: LD::is_r2 is set but never consumed, so r2 == r.
            assert_floats_close(actual_row, expected_row, output.with_suffix(".zld"))
        outputs[method] = output.with_suffix(".zld")
    assert outputs["r"].read_bytes() == outputs["r2"].read_bytes()
    return outputs["r"]


def solve_linear(matrix, vector):
    n = len(vector)
    augmented = [list(matrix[row]) + [vector[row]] for row in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        assert abs(augmented[pivot][column]) > 1e-12
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[-1] for row in augmented]


def residualize(values, design):
    columns = len(design[0])
    gram = [[sum(row[a] * row[b] for row in design) for b in range(columns)] for a in range(columns)]
    rhs = [sum(row[column] * value for row, value in zip(design, values)) for column in range(columns)]
    coefficients = solve_linear(gram, rhs)
    return [
        value - sum(coefficient * item for coefficient, item in zip(coefficients, row))
        for value, row in zip(values, design)
    ]


def write_lr_fixture(work):
    work = Path(work)
    markers = synthetic_markers(sample_count=12, marker_count=15)
    prefix = work / "input" / "lr"
    write_plink(prefix, markers)
    discrete = ["A" if index % 3 == 0 else "B" if index % 3 == 1 else "C" for index in range(12)]
    ranked = ["L" if index < 4 else "M" if index < 8 else "H" for index in range(12)]
    phenotypes = [
        0.4 * markers[1][index] - 0.2 * markers[2][index]
        + 0.3 * (index % 3) + 0.1 * index + math.sin(index)
        for index in range(12)
    ]
    pheno = work / "pheno.txt"
    pheno.write_text(
        "".join(f"F1 I{index} {value:.17g}\n" for index, value in enumerate(phenotypes)), encoding="utf-8"
    )
    covar = work / "covar.txt"
    covar.write_text(
        "".join(f"F1 I{index} {discrete[index]}\n" for index in range(12)), encoding="utf-8"
    )
    rcovar = work / "rcovar.txt"
    rcovar.write_text(
        "".join(f"F1 I{index} {ranked[index]}\n" for index in range(12)), encoding="utf-8"
    )
    design = []
    for index in range(12):
        # Any full-rank contrasts span the same residualizing subspace.
        design.append([
            1.0,
            float(discrete[index] == "B"), float(discrete[index] == "C"),
            float(ranked[index] == "M"), float(ranked[index] == "H"),
        ])
    return prefix, pheno, covar, rcovar, markers, phenotypes, design


def run_fastgwa_lr(gcta, work, threads):
    prefix, pheno, covar, rcovar, markers, phenotypes, design = write_lr_fixture(work)
    output = run_gcta(
        gcta,
        work,
        "lr",
        [
            "--bfile", prefix,
            "--fastGWA-lr",
            "--pheno", pheno,
            "--covar", covar,
            "--rcovar", rcovar,
            "--covar-maxlevel", 5,
            "--ge", 0.4, 0.6,
            "--verbose",
            "--seed", 12345,
        ],
        threads=threads,
    )
    assert_v2(output, "Performing fastGWA linear regression analysis")
    assert "Peak memory:" in log_text(output)
    assert "2 covar and 1 rcovar" not in log_text(output)
    rows = read_table(Path(f"{output}.fastGWA"))
    residual_y = residualize(phenotypes, design)
    expected_rows = [index for index, marker in enumerate(markers) if min(marker) != max(marker)]
    assert [int(row[1][2:]) for row in rows[1:]] == expected_rows
    for row, index in zip(rows[1:], expected_rows):
        residual_x = residualize([float(value) for value in markers[index]], design)
        expected_beta = sum(x * y for x, y in zip(residual_x, residual_y)) / sum(x * x for x in residual_x)
        assert math.isclose(float(row[7]), expected_beta, rel_tol=3e-5, abs_tol=2e-6), row
    return Path(f"{output}.fastGWA")


def make_family_fixture(work, sample_count=32, marker_count=131):
    markers, chromosomes = family_markers(sample_count=sample_count, marker_count=marker_count)
    prefix = Path(work) / "input" / "family"
    write_plink(prefix, markers, chromosomes)
    pheno, qcovar = write_family_phenotypes(prefix, markers)
    return prefix, pheno, qcovar, markers


def make_sparse(gcta, work, prefix, threads):
    return run_gcta(
        gcta, work, "sparse", ["--bfile", prefix, "--make-grm", "--sparse-cutoff", 0.05], threads=threads
    )


def run_fastgwa_exact(gcta, work, threads):
    work = Path(work)
    prefix, pheno, qcovar, markers = make_family_fixture(work)
    sparse = make_sparse(gcta, work, prefix, threads)
    common = [
        "--bfile", prefix,
        "--grm-sparse", sparse,
        "--fastGWA-mlm-exact",
        "--pheno", pheno,
        "--qcovar", qcovar,
        "--seed", 12345,
    ]
    direct = run_gcta(
        gcta,
        work,
        "direct",
        [*common, "--ge", 0.4, 0.6, "--ldlt"],
        threads=threads,
    )
    assert_v2(direct, "Inverting the variance-covariance matrix by ldlt")
    rows = read_table(Path(f"{direct}.fastGWA"))
    assert len(rows) == len(markers)  # Header plus all non-monomorphic markers.
    assert all(math.isfinite(float(row[7])) for row in rows[1:])

    saved = run_gcta(
        gcta, work, "saved", [*common, "--ge", 0.4, 0.6, "--save-inv"], threads=threads
    )
    assert_v2(saved, "The inverse of V matrix has been saved")
    inverse_prefix = Path(f"{saved}.fastGWA")
    assert Path(f"{inverse_prefix}.grm.inv").stat().st_size > 0
    loaded = run_gcta(
        gcta,
        work,
        "loaded",
        [*common, "--ge", 0.4, 0.6, "--load-inv", inverse_prefix],
        threads=threads,
    )
    assert_v2(loaded, "Inverse of V matrix loaded")
    assert_tables_close(Path(f"{direct}.fastGWA"), Path(f"{loaded}.fastGWA"))

    model = run_gcta(
        gcta,
        work,
        "reml",
        [
            *common,
            "--reml-detail",
            "--h2-limit", 1.2,
            "--gwa-no-constrain",
            "--no-marker",  # Currently accepted and listed, but not consumed by FastFAM.
            "--model-only",
        ],
        threads=threads,
    )
    assert_v2(model, "fastGWA-REML converged")
    assert Path(f"{model}.fastGWA.reml").stat().st_size > 0
    assert Path(f"{model}.fastGWA.mdl.bin").read_bytes()[:4] == b"fGWA"
    return Path(f"{direct}.fastGWA")


def run_fastgwa_grammar(gcta, work, threads):
    work = Path(work)
    prefix, pheno, qcovar, markers = make_family_fixture(work, sample_count=128, marker_count=2049)
    sparse = make_sparse(gcta, work, prefix, threads)
    output = run_gcta(
        gcta,
        work,
        "grammar",
        [
            "--bfile", prefix,
            "--grm-sparse", sparse,
            "--fastGWA-mlm",
            "--pheno", pheno,
            "--qcovar", qcovar,
            "--ge", 0.5, 0.5,
            "--save-fastGWA-mlm-residual",
            "--c-inf",
            "--c-inf-no-filter",
            "--seed", 12345,
        ],
        threads=threads,
    )
    assert_v2(output, "Performing fastGWA mixed model association analysis")
    assert "Mean GRAMMAR-Gamma value =" in log_text(output)
    gamma = float(Path(f"{output}.fastGWA.gamma").read_text(encoding="utf-8"))
    assert math.isfinite(gamma) and gamma > 0.0
    residual = read_table(Path(f"{output}.fastGWA.residual"))
    assert len(residual) == len(markers[0])
    assert all(math.isfinite(float(row[-1])) for row in residual)
    cinf = read_table(Path(f"{output}.fastGWA.cinf"))
    assert cinf[0] == ["CHR", "SNP", "POS", "A1", "A2", "cinf", "chisq"]
    assert len(cinf) >= 101
    assert len(read_table(Path(f"{output}.fastGWA"))) == len(markers)
    return Path(f"{output}.fastGWA")


def run_fastgwa_binary(gcta, work, threads):
    work = Path(work)
    prefix, pheno, qcovar, markers = make_family_fixture(work, sample_count=64, marker_count=2049)
    sparse = make_sparse(gcta, work, prefix, threads)
    output = run_gcta(
        gcta,
        work,
        "binary",
        [
            "--bfile", prefix,
            "--grm-sparse", sparse,
            "--fastGWA-mlm-binary",
            "--pheno", pheno,
            "--qcovar", qcovar,
            "--mpheno", 2,
            "--cv-threshold", 0.2,
            "--tao-start", 0.01,
            "--seed", 12345,
        ],
        threads=threads,
    )
    assert_v2(output, "Performing fastGWA generalized linear mixed model association analysis")
    assert "Initializing with tao:" in log_text(output)
    assert "reading 200 SNP" in log_text(output)
    rows = read_table(Path(f"{output}.fastGWA"))
    assert len(rows) == len(markers)
    assert rows[0][-1] == "CONVERGE"
    assert all(row[-1] == "1" for row in rows[1:])
    return Path(f"{output}.fastGWA")


def run_fastgwa_ge(gcta, work, threads):
    work = Path(work)
    prefix, pheno, qcovar, markers = make_family_fixture(work)
    sparse = make_sparse(gcta, work, prefix, threads)
    envir = work / "envir.txt"
    envir.write_text(
        "".join(f"F1 I{index} {((index * 7) % 13 - 6) / 3.0}\n" for index in range(len(markers[0]))),
        encoding="utf-8",
    )
    output = run_gcta(
        gcta,
        work,
        "ge",
        [
            "--bfile", prefix,
            "--grm-sparse", sparse,
            "--fastGWA-mlm",
            "--pheno", pheno,
            "--qcovar", qcovar,
            "--envir", envir,
            "--optimal-rho",
            "--noSandwich",
            "--grid-size", 3,
            "--est-vg", "HE",
            "--rel-only",
            "--force-gwa",
            "--seed", 12345,
        ],
        threads=threads,
    )
    assert_v2(output, "Performing fastGWA-GE mixed model association analysis using model-based variance estimator")
    assert "Use related pairs only." in log_text(output)
    assert "The optimal hyperparameter is" in log_text(output)
    rho = float(Path(f"{output}.fastGWA.optimalrho").read_text(encoding="utf-8"))
    assert rho in (0.0, 0.5, 1.0)
    rows = read_table(Path(f"{output}.fastGWA"))
    assert "BETA_INTER" in rows[0]
    assert len(rows) == len(markers)
    return Path(f"{output}.fastGWA")


def run_inversion_errors(gcta, work):
    work = Path(work)
    prefix, pheno, _, _ = make_family_fixture(work)
    sparse = make_sparse(gcta, work, prefix, 1)
    for flag in ("--cg", "--llt", "--pardiso", "--tcg", "--lscg", "--inv-t1"):
        expect_v2_error(
            gcta,
            work,
            flag[2:],
            [
                "--bfile", prefix,
                "--grm-sparse", sparse,
                "--fastGWA-mlm-exact",
                "--pheno", pheno,
                "--ge", 0.4, 0.6,
                flag,
            ],
            "unknown matrix inverse method.",
        )


def run_obsolete_errors(gcta, work):
    work = Path(work)
    prefix, pheno, _, _ = make_family_fixture(work)
    for flag in ("--fastGWA", "--fastGWA-gram"):
        expect_v2_error(
            gcta,
            work,
            flag[2:],
            ["--bfile", prefix, "--pheno", pheno, flag],
            "obsoleted flag. Try --fastGWA-mlm or --fastGWA-mlm-exact.",
        )


def run_burden(gcta, work, threads):
    work = Path(work)
    prefix, pheno, qcovar, markers = make_family_fixture(work)
    model = run_gcta(
        gcta,
        work,
        "model",
        [
            "--bfile", prefix,
            "--fastGWA-lr",
            "--pheno", pheno,
            "--qcovar", qcovar,
            "--model-only",
            "--seed", 12345,
        ],
        threads=threads,
    )
    genes = work / "genes.txt"
    genes.write_text("gene chr start end strand\nGENE1 1 1000 1100 +\n", encoding="utf-8")
    output = run_gcta(
        gcta,
        work,
        "burden",
        [
            "--bfile", prefix,
            "--load-model", Path(f"{model}.fastGWA"),
            "--set-list", genes,
            "--burden",
        ],
        threads=threads,
    )
    assert_v2(output, "Processing 1 genes by burden test")
    rows = read_table(Path(f"{output}.fastGWA"))
    assert rows[0] == ["GENE", "START", "END", "VAR_N", "MAF_MEAN", "BETA", "SE", "P"]
    assert len(rows) == 2 and rows[1][0] == "GENE1"
    assert int(rows[1][3]) > 0
    assert all(math.isfinite(float(value)) for value in rows[1][4:])
    return Path(f"{output}.fastGWA")


def run_gsmr(gcta, work, threads):
    prefix, exposures, outcomes = make_gsmr_fixture(work)
    output = run_gcta(
        gcta,
        work,
        "gsmr",
        [
            "--bfile", prefix,
            "--gsmr-file", exposures, outcomes,
            "--gsmr2-beta",
            "--diff-freq", 0.2,
            "--gwas-thresh", 5e-6,
            "--gsmr-snp", 5,
            "--gsmr-ld-fdr", 0.05,
            "--clump-kb", 10000,
            "--clump-r2", 0.9,
            "--heidi-thresh", 0,
        ],
        threads=threads,
    )
    assert_v2(output, "GSMR analyses completed.")
    rows = read_table(output.with_suffix(".gsmr"))
    # GSMR2 has a distinct estimator; this is a format/finiteness and thread
    # characterization, while the focused GSMR suite supplies GSMR1 oracles.
    assert len(rows) == 2 and rows[1][:2] == ["EXP", "OUT"]
    assert all(math.isfinite(float(value)) for value in rows[1][2:5])
    assert int(rows[1][5]) >= 5
    return output.with_suffix(".gsmr")


def run_gsmr_errors(gcta, work):
    work = Path(work)
    prefix, exposures, outcomes = make_gsmr_fixture(work)
    cases = (
        ("--gsmr-alg", 1, "--gsmr-alg has been superseded by --gsmr-direction."),
        (
            "--heidi-snp",
            10,
            "--heidi-snp is discontinued. Please use --gsmr-snp-min to specify minimum number of SNP instruments for the HEIDI-outlier analysis.",
        ),
        (
            "--clump-p1",
            5e-6,
            "--clump-p1 is discontinued. Please use --gwas-thresh to specify a p-value threshold for index SNPs.",
        ),
    )
    for flag, value, diagnostic in cases:
        expect_v2_error(
            gcta,
            work,
            flag[2:],
            ["--bfile", prefix, "--gsmr-file", exposures, outcomes, flag, value],
            diagnostic,
        )


def run_ecojo(gcta, work, threads):
    expression, summary, _ = write_ecojo_fixture(Path(work))
    output = run_gcta(
        gcta,
        work,
        "ecojo",
        ["--efile", expression, "--ecojo", summary, "--ecojo-slct", "--ecojo-p", 0.01],
        threads=threads,
    )
    assert_v2(output, "Options:")
    rows = read_table(Path(f"{output}.slct.ecojo"))
    assert len(rows) == 3
    assert [row[0] for row in rows[1:]] == ["P0", "P1"]
    assert all(math.isfinite(float(row[5])) for row in rows[1:])
    return Path(f"{output}.slct.ecojo")


def run_acat_wind(gcta, work):
    work = Path(work).resolve()
    genes, snps, _ = write_acat_inputs(work)
    output = run_gcta(
        gcta,
        work,
        "acat.csv",
        ["--acat", "--gene-list", genes, "--snp-list", snps, "--wind", 400, "--max-maf", 0.05],
    )
    rows = read_table(output)
    assert rows[1][:5] == ["1", "GENE1", "100", "200", "3"]
    assert 0.0 < float(rows[1][5]) < 1.0
    # ACAT returns before normal logger creation; its output and lack of the
    # legacy parser banner are the available dispatch evidence.
    assert "Accepted options:" not in console_text(output)


def run_allowlisted_noops(gcta, work):
    work = Path(work)
    markers = synthetic_markers(sample_count=8, marker_count=12)
    prefix = work / "input" / "noop"
    write_plink(prefix, markers)
    outputs = []
    for name, extra in (
        ("baseline", []),
        ("noblas", ["--noblas"]),
        ("num_vec", ["--num-vec", 10]),
        ("trace_exact", ["--trace-exact"]),
    ):
        output = run_gcta(
            gcta, work, name, ["--bfile", prefix, "--make-grm", *extra], threads=1
        )
        assert_v2(output, "The GRM computation is completed.")
        expected, counts = expected_grm(markers, 0)
        assert_floats_close(read_floats(Path(f"{output}.grm.bin")), expected, output)
        assert read_floats(Path(f"{output}.grm.N.bin")) == counts
        outputs.append(output)
    for output in outputs[1:]:
        assert_grm_equal(outputs[0], output)
    expect_v2_error(
        gcta,
        work,
        "bed_loader",
        [
            "--bim", prefix.with_suffix(".bim"),
            "--fam", prefix.with_suffix(".fam"),
            "--bed", prefix.with_suffix(".bed"),
            "--make-grm",
        ],
        "no genotype file is specified",
    )


def compare_thread_outputs(scenario, single, multi):
    if scenario in ("components", "pruning"):
        for first, second in zip(single, multi):
            if str(first).endswith(".frq"):
                assert first.read_bytes() == second.read_bytes()
            elif scenario == "components":
                assert_grm_equal(first, second)
            else:
                assert Path(f"{first}.grm.id").read_bytes() == Path(f"{second}.grm.id").read_bytes()
    elif scenario in ("multifile", "bgen_dosage"):
        assert_grm_equal(single, multi)
    elif scenario == "pgen_inputs":
        for first, second in zip(single, multi):
            assert_grm_equal(first, second)
    elif scenario == "grm_merge":
        assert Path(f"{single[0]}.grm.bin").read_bytes() == Path(f"{multi[0]}.grm.bin").read_bytes()
        assert_grm_equal(single[1], multi[1])
    elif scenario in ("filters", "ld", "partitions"):
        if isinstance(single, (tuple, list)):
            for first, second in zip(single, multi):
                assert first.read_bytes() == second.read_bytes()
        else:
            assert single.read_bytes() == multi.read_bytes()
    else:
        assert_tables_close(single, multi)


def main():
    parser = argparse.ArgumentParser(description="GCTA v2 supported-option behavioral tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    args = parser.parse_args()
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)

    error_only = {"inversion_errors", "obsolete_errors", "gsmr_errors", "allowlisted_noops"}
    if args.scenario == "acat_wind":
        run_acat_wind(args.gcta, work)
    elif args.scenario in error_only:
        runner = {
            "inversion_errors": run_inversion_errors,
            "obsolete_errors": run_obsolete_errors,
            "gsmr_errors": run_gsmr_errors,
            "allowlisted_noops": run_allowlisted_noops,
        }[args.scenario]
        runner(args.gcta, work)
    else:
        runners = {
            "components": run_components,
            "filters": run_filters,
            "pruning": run_pruning,
            "multifile": run_multifile,
            "pgen_inputs": run_pgen_inputs,
            "bgen_dosage": run_bgen_dosage,
            "grm_merge": run_grm_merge,
            "ld": run_ld,
            "fastgwa_lr": run_fastgwa_lr,
            "fastgwa_exact": run_fastgwa_exact,
            "fastgwa_grammar": run_fastgwa_grammar,
            "fastgwa_binary": run_fastgwa_binary,
            "fastgwa_ge": run_fastgwa_ge,
            "burden": run_burden,
            "gsmr": run_gsmr,
            "ecojo": run_ecojo,
        }
        if args.scenario == "partitions":
            single = (
                run_x_parts(args.gcta, work / "single_x", 1),
                run_dominance_parts(args.gcta, work / "single_d", 1),
            )
            multi = (
                run_x_parts(args.gcta, work / "multi_x", 4),
                run_dominance_parts(args.gcta, work / "multi_d", 4),
            )
            assert_grm_equal(single[0], multi[0])
            assert_floats_close(read_floats(single[1]), read_floats(multi[1]), multi[1])
        else:
            single = runners[args.scenario](args.gcta, work / "single", 1)
            multi = runners[args.scenario](args.gcta, work / "multi", 4)
            compare_thread_outputs(args.scenario, single, multi)
    print(f"Passed v2 option regression scenario: {args.scenario}")


if __name__ == "__main__":
    main()
