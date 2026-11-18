#!/usr/bin/env python3

import argparse
import gzip
import math
from pathlib import Path

from gcta_runner import assert_grm_equal, assert_tables_close, read_floats, read_table, run_gcta
from test_dosage import write_fixture as write_dosage_fixture
from test_errors import expect_error
from test_legacy_commands import (
    LD_MARKERS,
    assert_close,
    make_basic_fixture,
    make_family_grm,
    make_summary_fixture,
    path_with,
    run_legacy,
    write_snp_sets,
)
from test_simu import make_fixture as make_simu_fixture

SCENARIOS = (
    "filters",
    "grm_controls",
    "ld_controls",
    "reml_controls",
    "reml_inputs",
    "bivar_controls",
    "mlma_controls",
    "simu_controls",
    "dosage_controls",
    "misc_controls",
)


def frequency_ids(output):
    return [line.split()[0] for line in Path(output.with_suffix(".freq")).read_text().splitlines()]


def run_freq_case(gcta, work, name, prefix, extra, expected_ids, threads):
    output = run_legacy(
        gcta,
        work,
        name,
        ["--bfile", prefix, "--freq-v1", *extra],
        threads,
        "Calculating allele frequencies",
    )
    assert frequency_ids(output) == expected_ids
    return Path(output.with_suffix(".freq"))


def scenario_filters(gcta, work, threads):
    prefix, markers = make_basic_fixture(work)
    positions = [1000 + index * 5000 for index in range(len(markers))]
    prefix.with_suffix(".bim").write_text(
        "".join(f"1 rs{index} 0 {positions[index]} A G\n" for index in range(len(markers))),
        encoding="utf-8",
    )
    all_ids = [f"rs{index}" for index in range(len(markers))]
    cases = (
        ("extract-snp", ["--extract-snp", "rs3"], [all_ids[3]]),
        ("exclude-snp", ["--exclude-snp", "rs3"], [value for value in all_ids if value != "rs3"]),
        ("extract-region-snp", ["--extract-region-snp", "rs3", 2], [all_ids[3]]),
        ("extract-region-bp", ["--extract-region-bp", 1, positions[3], 2], [all_ids[3]]),
        ("exclude-region-snp", ["--exclude-region-snp", "rs3", 2], [value for value in all_ids if value != "rs3"]),
        ("exclude-region-bp", ["--exclude-region-bp", 1, positions[3], 2], [value for value in all_ids if value != "rs3"]),
    )
    results = []
    for name, extra, expected in cases:
        outputs = []
        for threads in (1, 4):
            outputs.append(run_freq_case(gcta, Path(work) / name / f"threads_{threads}", "freq", prefix, extra, expected, threads))
        assert outputs[0].read_bytes() == outputs[1].read_bytes()
        results.append(outputs[0])
    return results


def eigen_values(output):
    return [float(value) for value in Path(f"{output}.eigenval").read_text().splitlines()]


def assert_eigen_equal(first, second):
    values = eigen_values(first)
    reference = eigen_values(second)
    assert len(values) == len(reference)
    assert all(math.isclose(a, b, rel_tol=2e-4, abs_tol=2e-6) for a, b in zip(values, reference))


def scenario_grm_controls(gcta, work, threads):
    prefix, _ = make_basic_fixture(work)
    source = run_legacy(
        gcta, work, "source", ["--bfile", prefix, "--make-grm-v1"], threads, "GRM"
    )
    outputs = []

    dominance = run_legacy(
        gcta, work, "dominance", ["--bfile", prefix, "--make-grm-v1", "--dominance"], threads, "--dominance"
    )
    direct_dominance = run_legacy(
        gcta, work, "direct_dominance", ["--bfile", prefix, "--make-grm-d-v1"], threads, "GRM of"
    )
    assert Path(f"{dominance}.d.grm.bin").read_bytes() == Path(f"{direct_dominance}.d.grm.bin").read_bytes()
    assert Path(f"{dominance}.d.grm.N.bin").read_bytes() == Path(f"{direct_dominance}.d.grm.N.bin").read_bytes()
    outputs.append(Path(f"{dominance}.d.grm.bin"))

    adjusted = run_legacy(
        gcta, work, "adjusted", ["--grm", source, "--make-grm", "--grm-adj", 0.05], threads, "--grm-adj"
    )
    assert read_floats(Path(f"{adjusted}.grm.bin")) != read_floats(Path(f"{source}.grm.bin"))
    outputs.append(Path(f"{adjusted}.grm.bin"))

    high_ld = run_legacy(
        gcta, work, "high_ld", ["--bfile", prefix, "--make-grm-v1", "--rm-high-ld", 0.8], threads, "--rm-high-ld"
    )
    assert all(math.isfinite(value) for value in read_floats(Path(f"{high_ld}.grm.bin")))
    outputs.append(Path(f"{high_ld}.grm.bin"))

    cutoff = run_legacy(
        gcta, work, "cutoff", ["--grm", source, "--make-grm", "--grm-cutoff", 0.5, "--save-ram"], threads, "Pruning the GRM"
    )
    cutoff_alias = run_legacy(
        gcta, work, "cutoff_alias", ["--grm", source, "--make-grm", "--grm-cutoff-v1", 0.5, "--save-ram"], threads, "Pruning the GRM"
    )
    assert Path(f"{cutoff}.grm.id").read_bytes() == Path(f"{cutoff_alias}.grm.id").read_bytes()
    assert Path(f"{cutoff}.grm.bin").read_bytes() == Path(f"{cutoff_alias}.grm.bin").read_bytes()
    outputs.append(Path(f"{cutoff_alias}.grm.bin"))

    binary_pca = run_legacy(gcta, work, "binary_pca", ["--grm-bin", source, "--pca", 2], threads, "Performing principal component analysis")
    baseline_pca = run_legacy(gcta, work, "baseline_pca", ["--grm", source, "--pca", 2], threads, "Performing principal component analysis")
    assert_eigen_equal(binary_pca, baseline_pca)
    outputs.append(Path(f"{binary_pca}.eigenval"))

    compressed = run_legacy(
        gcta, work, "compressed", ["--bfile", prefix, "--make-grm-gz"], threads, "GRM"
    )
    gz_pca = run_legacy(
        gcta, work, "gz_pca", ["--grm-gz", compressed, "--pca", 2], threads, "Performing principal component analysis"
    )
    assert_eigen_equal(gz_pca, baseline_pca)
    outputs.append(Path(f"{gz_pca}.eigenval"))

    binary_list = Path(work) / "binary_grms.txt"
    binary_list.write_text(f"{source}\n", encoding="utf-8")
    mgrm_bin = run_legacy(gcta, work, "mgrm_bin", ["--mgrm-bin", binary_list, "--pca", 2], threads, "Performing principal component analysis")
    assert_eigen_equal(mgrm_bin, baseline_pca)
    outputs.append(Path(f"{mgrm_bin}.eigenval"))

    gz_list = Path(work) / "gz_grms.txt"
    gz_list.write_text(f"{compressed}\n", encoding="utf-8")
    mgrm_gz = run_legacy(gcta, work, "mgrm_gz", ["--mgrm-gz", gz_list, "--pca", 2], threads, "Performing principal component analysis")
    assert_eigen_equal(mgrm_gz, baseline_pca)
    outputs.append(Path(f"{mgrm_gz}.eigenval"))

    prune_default = run_legacy(
        gcta, work, "prune_default", ["--bfile", prefix, "--ld-pruning", 0.5, "--ld-wind", 1], threads, "Pruning SNPs for LD"
    )
    prune_alias = run_legacy(
        gcta,
        work,
        "prune_alias",
        ["--bfile", prefix, "--ld-pruning", 0.5, "--make-grm-wt-wind", 1],
        1,
        "Pruning SNPs for LD",
    )
    assert Path(f"{prune_default}.prune.in").read_bytes() == Path(f"{prune_alias}.prune.in").read_bytes()
    outputs.append(Path(f"{prune_alias}.prune.in"))
    return outputs


def scenario_ld_controls(gcta, work, threads):
    prefix = Path(work) / "input" / "ld"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    from plink import write_plink

    write_plink(prefix, LD_MARKERS)
    target = Path(work) / "target.txt"
    target.write_text("rs0\n", encoding="utf-8")
    outputs = []

    base = run_legacy(
        gcta, work, "base", ["--bfile", prefix, "--ld", target, "--ld-wind", 1000], threads, "Estimating LD structure"
    )
    step = run_legacy(
        gcta,
        work,
        "step",
        ["--bfile", prefix, "--ld", target, "--ld-wind", 1000, "--ld-step", 1],
        1,
        "Estimating LD structure",
    )
    assert Path(f"{step}.rsq.ld").exists() and Path(f"{base}.rsq.ld").exists()
    outputs.append(Path(f"{step}.rsq.ld"))

    include_target = run_legacy(
        gcta,
        work,
        "include_target",
        ["--bfile", prefix, "--ld", target, "--ld-wind", 1000, "--ld-i"],
        1,
        "Estimating LD structure",
    )
    assert "rs0" in Path(f"{include_target}.snp.ld").read_text().split()
    assert "rs0" not in Path(f"{base}.snp.ld").read_text().split()
    outputs.append(Path(f"{include_target}.snp.ld"))

    prune_default = run_legacy(
        gcta, work, "prune_default", ["--bfile", prefix, "--ld-pruning", 0.5, "--ld-wind", 1000], threads, "Pruning SNPs for LD"
    )
    prune_alias = run_legacy(
        gcta,
        work,
        "prune_alias",
        ["--bfile", prefix, "--ld-pruning", 0.5, "--ld-pruning-wind", 1000],
        1,
        "Pruning SNPs for LD",
    )
    assert Path(f"{prune_default}.prune.in").read_bytes() == Path(f"{prune_alias}.prune.in").read_bytes()
    outputs.append(Path(f"{prune_alias}.prune.in"))

    score = run_legacy(
        gcta, work, "score", ["--bfile", prefix, "--ld-score", "--ld-wind", 1000], threads, "Calculating LD score"
    )
    score_filtered = run_legacy(
        gcta,
        work,
        "score_filtered",
        ["--bfile", prefix, "--ld-score", "--ld-wind", 1000, "--ld-rsq-cutoff", 0.9],
        1,
        "Calculating LD score",
    )
    score_adjusted = run_legacy(
        gcta,
        work,
        "score_adjusted",
        ["--bfile", prefix, "--ld-score", "--ld-wind", 1000, "--ld-score-adj"],
        1,
        "Calculating LD score",
    )
    for path in (f"{score}.score.ld", f"{score_filtered}.score.ld", f"{score_adjusted}.score.ld"):
        rows = [row for row in read_table(path) if row]
        assert rows and all(math.isfinite(float(value)) for row in rows[1:] for value in row[3:7])
    outputs.extend((Path(f"{score_filtered}.score.ld"), Path(f"{score_adjusted}.score.ld")))

    return outputs


def scenario_reml_controls(gcta, work, threads):
    grm, pheno = make_family_grm(gcta, work, sample_count=64, marker_count=201)
    base = ["--grm", grm, "--pheno", pheno, "--reml", "--reml-priors", 0.5, "--reml-maxit", 100]
    variants = (
        ("force-converge", ["--reml-no-lrt", "--reml-force-converge"]),
        ("allow-no-converge", ["--reml-no-lrt", "--reml-allow-no-converge"]),
        ("amzvc", ["--reml-no-lrt", "--reml-amzvc"]),
        ("bending", ["--reml-no-lrt", "--reml-bending"]),
        ("bendV", ["--reml-no-lrt", "--reml-bendV"]),
        ("diag-one", ["--reml-no-lrt", "--reml-diag-one"]),
        ("diag-adj", ["--reml-no-lrt", "--reml-diagV-adj", 1]),
        ("diag-mul", ["--reml-no-lrt", "--reml-diag-mul", 0.05]),
        ("inv-mtd", ["--reml-no-lrt", "--reml-inv-mtd", 1]),
        ("pred-rand", ["--reml-no-lrt", "--reml-pred-rand"]),
        ("est-fix-varcov", ["--reml-est-fix-varcov"]),
        ("priors-var", ["--reml-no-lrt", "--reml-priors-var", 0.5, 0.5]),
        ("fixed-var", ["--reml-no-lrt", "--reml-fixed-var", 0.5, 0.5]),
        ("lrt", ["--reml-lrt", 1]),
    )
    outputs = []
    for name, flags in variants:
        case = Path(work) / name
        output = run_legacy(gcta, case, "reml", [*base, *flags], threads, "Analysis finished")
        hsq = Path(f"{output}.hsq")
        assert hsq.exists(), name
        values = {}
        for row in read_table(hsq):
            if len(row) >= 2:
                try:
                    values[row[0]] = float(row[1])
                except ValueError:
                    pass
        assert values["n"] == 64
        assert all(math.isfinite(values[key]) for key in ("V(G)", "V(e)", "Vp", "V(G)/Vp", "logL"))
        outputs.append(hsq)
    return outputs



def scenario_reml_inputs(gcta, work, threads):
    grm, pheno = make_family_grm(gcta, work, sample_count=64, marker_count=201)
    phenotypes = [line.split() for line in Path(pheno).read_text(encoding="utf-8").splitlines()]
    discrete = Path(work) / "gxe.txt"
    discrete.write_text(
        "".join(f"{row[0]} {row[1]} {'A' if index % 2 else 'B'}\n" for index, row in enumerate(phenotypes)),
        encoding="utf-8",
    )
    quantitative = Path(work) / "gxqe.txt"
    quantitative.write_text(
        "".join(f"{row[0]} {row[1]} {(index % 7 - 3) / 3.0}\n" for index, row in enumerate(phenotypes)),
        encoding="utf-8",
    )
    weights = Path(work) / "weights.txt"
    weights.write_text(
        "".join(f"{row[0]} {row[1]} 1.0\n" for row in phenotypes), encoding="utf-8"
    )
    variants = (
        ("prevalence", ["--reml-priors", 0.5, "--mpheno", 2, "--prevalence", 0.1]),
        ("cvblup", ["--reml-priors", 0.5, "--cvblup"]),
        ("gxe", ["--reml-priors", 0.5, 0.5, "--gxe", discrete]),
        ("gxqe", ["--reml-priors", 0.5, 0.5, "--gxqe", quantitative]),
        ("res-diag", ["--reml-priors", 0.5, "--reml-res-diag", weights]),
    )
    base = ["--grm", grm, "--pheno", pheno, "--reml", "--reml-maxit", 100, "--reml-no-lrt"]
    outputs = []
    for name, flags in variants:
        case = Path(work) / name
        output = run_legacy(gcta, case, "reml", [*base, *flags], threads, "Analysis finished")
        hsq = Path(f"{output}.hsq")
        assert hsq.exists() and hsq.stat().st_size > 0, name
        assert any(math.isfinite(float(row[1])) for row in read_table(hsq) if len(row) > 1 and row[0] in ("V(G)", "V(G)/Vp"))
        outputs.append(hsq)
    return outputs


def scenario_bivar_controls(gcta, work, threads):
    from test_bivar import make_fixture as make_bivar_fixture

    grm, pheno = make_bivar_fixture(work, gcta)
    common = ["--grm", grm, "--pheno", pheno, "--reml-bivar", "--reml-maxit", 100]
    variants = (
        ("baseline", ["--reml-no-lrt"]),
        ("no-constrain", ["--reml-no-lrt", "--reml-bivar-no-constrain"]),
        ("nocove", ["--reml-no-lrt", "--reml-bivar-nocove"]),
        ("prevalence", ["--reml-no-lrt", "--reml-bivar-prevalence", 0.1, 0.1]),
        ("lrt-rg", ["--reml-bivar-lrt-rg", 0]),
    )
    outputs = []
    for name, flags in variants:
        output = run_legacy(gcta, Path(work) / name, "bivar", [*common, *flags], threads, "Log-likelihood ratio converged")
        hsq = Path(f"{output}.hsq")
        values = {}
        for row in read_table(hsq):
            if len(row) >= 2:
                try:
                    values[row[0]] = float(row[1])
                except ValueError:
                    pass
        assert values["n"] == 256
        assert -1.5 <= values["rG"] <= 1.5
        assert math.isfinite(values["logL"])
        outputs.append(hsq)
    return outputs


def scenario_mlma_controls(gcta, work, threads):
    grm, pheno = make_family_grm(gcta, work, sample_count=64, marker_count=201)
    qcovar = Path(grm).with_suffix(".qcovar")
    if not qcovar.exists():
        qcovar = Path(pheno).with_suffix(".qcovar")
    source_prefix = Path(pheno).with_suffix("")
    extract = Path(work) / "subset.snplist"
    extract.write_text(
        "".join(f"rs{index}\n" for index in range(0, 201, 2)), encoding="utf-8"
    )
    subset_grm = run_legacy(
        gcta, work, "subset_grm", ["--bfile", source_prefix, "--extract", extract, "--make-grm-v1"], threads, "GRM"
    )
    base = ["--bfile", source_prefix, "--grm", grm, "--pheno", pheno, "--mlma", "--reml-priors", 0.5, "--reml-maxit", 100]
    variants = (
        ("baseline", []),
        ("no-adj", ["--mlma-no-adj-covar"]),
        ("no-preadj", ["--mlma-no-preadj-covar"]),
        ("subtract", ["--mlma-subtract-grm", subset_grm]),
    )
    outputs = []
    for name, flags in variants:
        arguments = [*base, *flags]
        output = run_legacy(gcta, Path(work) / name, "mlma", arguments, threads, "Running association tests")
        rows = [row for row in read_table(Path(f"{output}.mlma")) if row]
        assert len(rows) == 202
        assert all(math.isfinite(float(row[8])) for row in rows[1:] if 0.0 < float(row[5]) < 1.0)
        outputs.append(Path(f"{output}.mlma"))
    assert [row[1] for row in read_table(outputs[1])[1:]] == [row[1] for row in read_table(outputs[0])[1:]]
    return outputs


def scenario_simu_controls(gcta, work, threads):
    prefix, qtl = make_simu_fixture(work)
    output = run_legacy(
        gcta,
        work,
        "simu",
        [
            "--bfile", prefix,
            "--simu-qt",
            "--simu-rep", 1,
            "--simu-hsq", 0.5,
            "--simu-causal-loci", qtl,
            "--simu-seed", 12345,
            "--simu-eff-mod", 1,
            "--simu-embayesb",
            "--simu-ouput-causal",
        ],
        threads,
        "Simulating",
    )
    par_rows = read_table(Path(f"{output}.par"))
    phen_rows = read_table(Path(f"{output}.emb"))
    assert par_rows[0] == ["QTL", "RefAllele", "Frequency", "Effect"]
    assert [row[0] for row in par_rows[1:]] == ["rs1", "rs2"]
    assert all(math.isfinite(float(row[3])) for row in par_rows[1:])
    assert len(phen_rows) >= 128
    assert all(math.isfinite(float(value)) for row in phen_rows for value in row[1:])
    return Path(f"{output}.par")


def scenario_dosage_controls(gcta, work, threads):
    dose, info = write_dosage_fixture(work, "mach")
    baseline = run_legacy(
        gcta, work, "mach", ["--dosage-mach", dose, info, "--make-grm"], threads, "GRM"
    )
    compressed = Path(f"{dose}.gz")
    with gzip.open(compressed, "wt", encoding="utf-8") as stream:
        stream.write(Path(dose).read_text(encoding="utf-8"))
    gz = run_legacy(
        gcta, work, "mach_gz", ["--dosage-mach-gz", compressed, info, "--make-grm"], threads, "GRM"
    )
    assert_grm_equal(baseline, gz)

    info_rows = Path(info).read_text(encoding="utf-8").splitlines()
    fields = info_rows[1].split()
    fields[-1] = "0.5"
    info_rows[1] = " ".join(fields)
    filtered_info = Path(work) / "filtered.info"
    filtered_info.write_text("\n".join(info_rows) + "\n", encoding="utf-8")
    filtered = run_legacy(
        gcta,
        work,
        "filtered",
        ["--dosage-mach", dose, filtered_info, "--imput-rsq", 0.8, "--make-grm"],
        threads,
        "GRM",
    )
    assert "After filtering for imputation Rsq > 0.8, there are 3 SNPs" in Path(f"{filtered}.log").read_text(encoding="utf-8")

    updated_info = Path(work) / "updated_rsq.txt"
    updated_info.write_text("".join(f"rs{index} 0.95\n" for index in range(4)), encoding="utf-8")
    updated = run_legacy(
        gcta,
        work,
        "updated",
        ["--dosage-mach", dose, filtered_info, "--update-imput-rsq", updated_info, "--make-grm"],
        threads,
        "GRM",
    )
    assert all(math.isfinite(value) for value in read_floats(Path(f"{updated}.grm.bin")))
    return Path(f"{gz}.grm.bin"), Path(f"{filtered}.grm.bin"), Path(f"{updated}.grm.bin")


def scenario_misc_controls(gcta, work, threads):
    prefix, _ = make_basic_fixture(work)
    expect_error(
        gcta,
        work,
        "bfile2_without_bfile",
        ["--bfile2", prefix, "--freq-v1", "--out", Path(work) / "result"],
    )
    assert "the option --bfile2 should always go with the option --bfile" in Path(
        work / "bfile2_without_bfile.console"
    ).read_text(encoding="utf-8")
    freq = run_legacy(gcta, work, "freq", ["--bfile", prefix, "--freq-v1"], threads, "Calculating allele frequencies")
    ssq = run_legacy(
        gcta, work, "ssq", ["--bfile", prefix, "--freq-v1", "--ssq"], threads, "Calculating allele frequencies"
    )
    assert Path(f"{ssq}.freq").exists()
    assert Path(f"{ssq}.freq").read_text() == Path(f"{freq}.freq").read_text()
    return Path(f"{ssq}.freq")

def main():
    parser = argparse.ArgumentParser(description="Legacy modifier regression tests")
    parser.add_argument("--gcta", required=True)
    parser.add_argument("--data")
    parser.add_argument("--work", required=True)
    parser.add_argument("--scenario", required=True, choices=SCENARIOS)
    args = parser.parse_args()
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    scenarios = {
        "filters": scenario_filters,
        "grm_controls": scenario_grm_controls,
        "ld_controls": scenario_ld_controls,
        "reml_controls": scenario_reml_controls,
        "reml_inputs": scenario_reml_inputs,
        "bivar_controls": scenario_bivar_controls,
        "mlma_controls": scenario_mlma_controls,
        "simu_controls": scenario_simu_controls,
        "dosage_controls": scenario_dosage_controls,
        "misc_controls": scenario_misc_controls,
    }
    if args.scenario not in scenarios:
        raise NotImplementedError(args.scenario)
    single = scenarios[args.scenario](args.gcta, work / "single", 1)
    multi = scenarios[args.scenario](args.gcta, work / "multi", 4)
    if isinstance(single, Path):
        single = (single,)
    if isinstance(multi, Path):
        multi = (multi,)
    for first, second in zip(single, multi):
        assert Path(first).read_bytes() == Path(second).read_bytes(), (first, second)
    print(f"Passed legacy modifier scenario: {args.scenario}")


if __name__ == "__main__":
    main()
