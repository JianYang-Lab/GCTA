#!/usr/bin/env python3
"""
MCP Server for GCTA (Genome-wide Complex Trait Analysis).

Exposes all GCTA functionality through MCP tools so that AI assistants
like Claude, Codex, and opencode can run GCTA analyses via natural language.

Usage:
    python server.py

Configuration via environment variables:
    GCTA_BINARY_PATH  Path to the GCTA binary. Defaults to
                      build/Release/gcta64 in the source tree, or ../gcta
                      relative to this script in a release package.
    GCTA_WORK_DIR     Working directory for GCTA data files
                      (default: source or release package root)
    GCTA_TIMEOUT      Execution timeout in seconds (default: 3600)
"""

import os
import sys
import json
import subprocess
import glob

from mcp.server.fastmcp import FastMCP


# ============================================================================
# Configuration
# ============================================================================

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PKG_ROOT = os.path.dirname(_SCRIPT_DIR)
_IS_SOURCE_TREE = os.path.basename(_PKG_ROOT) == ".opencode"
_SOURCE_ROOT = os.path.dirname(_PKG_ROOT) if _IS_SOURCE_TREE else _PKG_ROOT
_DEFAULT_BINARY = (
    os.path.join(_SOURCE_ROOT, "build", "Release", "gcta64")
    if _IS_SOURCE_TREE
    else os.path.join(_PKG_ROOT, "gcta")
)

GCTA_BINARY_PATH = os.environ.get("GCTA_BINARY_PATH", _DEFAULT_BINARY)
GCTA_WORK_DIR = os.environ.get("GCTA_WORK_DIR", _SOURCE_ROOT)
GCTA_TIMEOUT = int(os.environ.get("GCTA_TIMEOUT", "3600"))


# ============================================================================
# Execution wrapper
# ============================================================================

class GCTAExecutor:
    """Handles execution of GCTA binary via subprocess."""

    def __init__(self, binary_path, work_dir, timeout):
        self.binary_path = binary_path
        self.work_dir = work_dir
        self.timeout = timeout

    def execute(self, args: str, work_dir: str = None) -> dict:
        """Execute GCTA with the given command-line arguments."""
        wd = os.path.abspath(work_dir or self.work_dir)
        cmd = [self.binary_path] + args.split()
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=wd,
                timeout=self.timeout,
            )
            return {
                "exit_code": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": " ".join(cmd),
                "work_dir": wd,
            }
        except subprocess.TimeoutExpired:
            return {"error": f"GCTA execution timed out after {self.timeout}s"}
        except FileNotFoundError:
            return {"error": f"GCTA binary not found at {self.binary_path}"}
        except Exception as e:
            return {"error": str(e)}


# ============================================================================
# Helper functions
# ============================================================================

def _parse_out_prefix(args: str) -> str:
    """Parse the --out value from args; return default if absent."""
    tokens = args.split()
    for i, t in enumerate(tokens):
        if t == "--out" and i + 1 < len(tokens):
            return tokens[i + 1]
    return "gcta"


def _read_log_file(work_dir: str, out_prefix: str) -> str:
    """Read the GCTA .log file if it exists."""
    log_path = os.path.join(work_dir, out_prefix + ".log")
    if os.path.exists(log_path):
        try:
            with open(log_path, "r", errors="replace") as f:
                return f.read()
        except Exception:
            return ""
    return ""


def _list_output_files(work_dir: str, out_prefix: str) -> list:
    """List output files matching the out prefix."""
    files = []
    pattern = os.path.join(work_dir, out_prefix + "*")
    for f in sorted(glob.glob(pattern)):
        files.append({
            "file": os.path.relpath(f, work_dir),
            "size_bytes": os.path.getsize(f),
        })
    return files


def _format_result(exec_result: dict, work_dir: str,
                   out_prefix: str) -> str:
    """Format the execution result as a JSON string."""
    if "error" in exec_result:
        return json.dumps({"error": exec_result["error"]}, indent=2)

    log_content = _read_log_file(work_dir, out_prefix)
    output_files = _list_output_files(work_dir, out_prefix)

    result = {
        "exit_code": exec_result["exit_code"],
        "success": exec_result["exit_code"] == 0,
        "command": exec_result["command"],
        "work_dir": work_dir,
        "out_prefix": out_prefix,
        "stdout": exec_result["stdout"][-5000:] if exec_result["stdout"] else "",
        "stderr": exec_result["stderr"][-5000:] if exec_result["stderr"] else "",
        "log_content": log_content[-8000:] if log_content else "",
        "output_files": output_files,
    }
    return json.dumps(result, indent=2, ensure_ascii=False)


def _build_args(base: list, **kwargs) -> str:
    """Build a GCTA argument string from a base list and keyword arguments.

    ``base`` is a list of required flags + values (e.g.
    ``["--bfile", bfile, "--make-grm", "--out", out]``).

    Keyword arguments whose value is:
      - ``None`` or ``""``  -> skipped
      - ``False``           -> skipped (boolean flags are added only when True)
      - ``True``            -> added as a bare flag (``--flag``)
      - any other value     -> added as ``--flag value``
    """
    parts = list(base)
    for key, val in kwargs.items():
        if val is None or val == "" or val is False:
            continue
        flag = "--" + key.replace("_", "-")
        if val is True:
            parts.append(flag)
        else:
            parts.extend([flag, str(val)])
    return " ".join(parts)


# ============================================================================
# Create MCP server
# ============================================================================

mcp = FastMCP("gcta-mcp")

_executor = GCTAExecutor(
    binary_path=GCTA_BINARY_PATH,
    work_dir=GCTA_WORK_DIR,
    timeout=GCTA_TIMEOUT,
)


# ============================================================================
# Core tools
# ============================================================================

@mcp.tool()
def gcta_run(args: str, work_dir: str = "") -> str:
    """Execute GCTA with arbitrary command-line arguments.

    This is the most flexible tool — it gives full access to every GCTA
    option.  Pass everything that would appear after ``gcta64`` on the
    command line.

    Examples:
        args = "--bfile test --make-grm --out test_grm"
        args = "--grm test_grm --reml --pheno test.phen --out test_reml"
        args = "--bfile test --mlma-loco --pheno test.phen --out test_mlma"

    Args:
        args: Command-line arguments (everything after the binary name).
        work_dir: Working directory for GCTA.  Defaults to GCTA_WORK_DIR.
                  All input/output file paths are relative to this directory.

    Returns:
        JSON with exit_code, stdout, stderr, log_content, output_files.
    """
    wd = work_dir or GCTA_WORK_DIR
    out_prefix = _parse_out_prefix(args)
    result = _executor.execute(args, wd)
    return _format_result(result, os.path.abspath(wd), out_prefix)


@mcp.tool()
def gcta_help(topic: str = "all") -> str:
    """Get comprehensive documentation about GCTA analyses and options.

    Args:
        topic: One of:
            "all"             - Overview of all analyses (default)
            "data"            - Data input, management, and filtering
            "grm"             - Genetic Relationship Matrix (GRM) construction
            "reml"            - REML / GREML variance component estimation
            "bivariate_reml"  - Bivariate REML analysis
            "hereg"           - Haseman-Elston regression
            "pca"             - Principal component analysis
            "mlma"            - Mixed linear model association (MLMA / MLMA-LOCO)
            "cojo"            - Conditional and joint analysis (COJO)
            "gsmr"            - GSMR Mendelian randomization
            "mtcojo"          - Multi-trait COJO
            "fastgwa"         - fastGWA genome-wide association
            "fastbat"         - fastBAT / mBAT gene-based test
            "simulation"      - Phenotype simulation
            "ld"              - LD pruning, LD score regression
            "fst"             - Fst population differentiation
            "expression"      - Expression data analysis (eRFile, ecojo, make-erm)
            "acat"            - ACAT gene-based test
            "options"         - Full list of all command-line options

    Returns:
        Documentation string for the requested topic.
    """
    docs = _get_help_docs()
    topic = topic.lower().strip()
    if topic in docs:
        return docs[topic]
    elif topic == "all":
        return docs["all"]
    else:
        available = ", ".join(sorted(docs.keys()))
        return f"Unknown topic '{topic}'. Available topics: {available}"


@mcp.tool()
def gcta_check() -> str:
    """Check if GCTA is properly configured and can be executed.

    Verifies:
      - GCTA binary exists at the configured path
      - Working directory exists and is writable
      - Test data files are present

    Returns:
        JSON with configuration details and check results.
    """
    info = {
        "binary_path": GCTA_BINARY_PATH,
        "binary_exists": os.path.exists(GCTA_BINARY_PATH),
        "work_dir": GCTA_WORK_DIR,
        "work_dir_exists": os.path.exists(GCTA_WORK_DIR),
        "work_dir_writable": os.access(GCTA_WORK_DIR, os.W_OK) if os.path.exists(GCTA_WORK_DIR) else False,
        "timeout_seconds": GCTA_TIMEOUT,
    }

    checks = []

    if not info["binary_exists"]:
        checks.append("FAIL: GCTA binary not found. Set GCTA_BINARY_PATH env var.")
    else:
        checks.append("OK: GCTA binary found.")

    if not info["work_dir_exists"]:
        checks.append("FAIL: Working directory not found. Set GCTA_WORK_DIR env var.")
    else:
        checks.append("OK: Working directory exists.")
        if info["work_dir_writable"]:
            checks.append("OK: Working directory is writable.")
        else:
            checks.append("WARN: Working directory is not writable.")

    test_data_dir = GCTA_WORK_DIR
    test_files = ["test.bed", "test.bim", "test.fam", "test.phen"]
    missing = [f for f in test_files if not os.path.exists(os.path.join(test_data_dir, f))]
    if not missing:
        checks.append("OK: Test data files found in working directory.")
    else:
        checks.append(f"WARN: Missing test data files: {', '.join(missing)}")

    info["checks"] = checks
    return json.dumps(info, indent=2, ensure_ascii=False)


# ============================================================================
# GRM construction tools
# ============================================================================

@mcp.tool()
def gcta_make_grm(
    bfile: str,
    out: str = "gcta_grm",
    maf: float = 0.0,
    thread_num: int = 1,
    autosome_num: int = 22,
    chr_filter: int = 0,
    autosome: bool = False,
    make_grm_alg: int = 0,
    grm_adj: float = -1.0,
    grm_cutoff: float = -1.0,
    dosage_compen: int = -1,
    dominance: bool = False,
    make_grm_xchr: bool = False,
    make_grm_inbred: bool = False,
    save_ram: bool = False,
    extract: str = "",
    exclude: str = "",
    work_dir: str = "",
) -> str:
    """Construct a Genetic Relationship Matrix (GRM) from PLINK bed data.

    The GRM estimates pairwise genomic relatedness between individuals
    based on genome-wide SNPs.  This is the foundational step for most
    downstream GCTA analyses (REML, PCA, MLMA, etc.).

    Args:
        bfile: PLINK binary file prefix (e.g. "test" for test.bed/.bim/.fam).
        out: Output file prefix.  Produces {out}.grm.bin, {out}.grm.N.bin,
             {out}.grm.id.
        maf: Minor allele frequency filter (0 = no filter, typical: 0.01).
        thread_num: Number of CPU threads.
        autosome_num: Number of autosomes (default 22 for human).
        chr_filter: Restrict to a specific chromosome (0 = no restriction).
        autosome: If True, restrict to autosomes 1..autosome_num.
        make_grm_alg: Algorithm: 0 = VanRaden (default), 1 = alternative.
        grm_adj: GRM adjustment factor (-1 = not used).
        grm_cutoff: Remove pairs with GRM > cutoff (-1 = no cutoff).
        dosage_compen: Dosage compensation: 0 = none, 1 = X chromosome.
        dominance: If True, calculate dominance GRM.
        make_grm_xchr: If True, calculate X-chromosome GRM.
        make_grm_inbred: If True, calculate inbreeding GRM.
        save_ram: If True, use memory-saving mode.
        extract: File of SNP IDs to include.
        exclude: File of SNP IDs to exclude.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--autosome-num", str(autosome_num),
            "--thread-num", str(thread_num)]

    if dominance:
        base.append("--make-grm-d")
    elif make_grm_xchr:
        base.append("--make-grm-xchr")
    elif make_grm_inbred:
        base.append("--make-grm-inbred")
    else:
        base.append("--make-grm")

    base.extend(["--out", out])

    args = _build_args(base, maf=maf, chr=chr_filter if chr_filter > 0 else None,
                       autosome=autosome, make_grm_alg=make_grm_alg,
                       grm_adj=grm_adj if grm_adj >= 0 else None,
                       grm_cutoff=grm_cutoff if grm_cutoff > 0 else None,
                       dc=dosage_compen if dosage_compen >= 0 else None,
                       save_ram=save_ram, extract=extract or None,
                       exclude=exclude or None)
    return gcta_run(args, work_dir)


# ============================================================================
# REML / GREML tools
# ============================================================================

@mcp.tool()
def gcta_reml(
    grm: str,
    pheno: str,
    out: str = "gcta_reml",
    qcovar: str = "",
    covar: str = "",
    mpheno: int = 0,
    thread_num: int = 1,
    reml_priors: str = "",
    reml_priors_var: str = "",
    reml_maxit: int = 100,
    prevalence: float = -1.0,
    reml_no_constrain: bool = False,
    reml_no_lrt: bool = False,
    reml_pred_rand: bool = False,
    reml_est_fix: bool = False,
    grm_cutoff: float = -1.0,
    gxe: str = "",
    work_dir: str = "",
) -> str:
    """Run REML analysis to estimate variance components (GREML).

    Estimates the proportion of phenotypic variance explained by all
    genome-wide SNPs (SNP-based heritability, h2_SNP) using Restricted
    Maximum Likelihood (REML).

    Args:
        grm: GRM file prefix (from gcta_make_grm).
        pheno: Phenotype file path.
        out: Output file prefix.  Produces {out}.hsq.
        qcovar: Quantitative covariate file.
        covar: Categorical covariate file.
        mpheno: Phenotype column number (0 = first column).
        thread_num: Number of CPU threads.
        reml_priors: Prior values of variance explained (space-separated).
        reml_priors_var: Prior values of variance components (space-separated).
        reml_maxit: Maximum REML iterations (default 100).
        prevalence: Disease prevalence for liability conversion (-1 = not used).
        reml_no_constrain: If True, do not constrain variance components.
        reml_no_lrt: If True, skip likelihood ratio test.
        reml_pred_rand: If True, predict random effects (BLUP).
        reml_est_fix: If True, estimate fixed effects.
        grm_cutoff: GRM cutoff for related pairs (-1 = no cutoff).
        gxe: GxE interaction file.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--grm", grm, "--pheno", pheno,
            "--thread-num", str(thread_num), "--reml",
            "--reml-maxit", str(reml_maxit), "--out", out]

    args = _build_args(base, qcovar=qcovar or None, covar=covar or None,
                       mpheno=mpheno or None,
                       reml_priors=reml_priors or None,
                       reml_priors_var=reml_priors_var or None,
                       prevalence=prevalence if prevalence > 0 else None,
                       reml_no_constrain=reml_no_constrain,
                       reml_no_lrt=reml_no_lrt,
                       reml_pred_rand=reml_pred_rand,
                       reml_est_fix=reml_est_fix,
                       grm_cutoff=grm_cutoff if grm_cutoff > 0 else None,
                       gxe=gxe or None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_bivariate_reml(
    grm: str,
    pheno: str,
    out: str = "gcta_bireml",
    qcovar: str = "",
    covar: str = "",
    mphen: int = 1,
    mphen2: int = 2,
    thread_num: int = 1,
    reml_maxit: int = 100,
    reml_bivar_prevalence: str = "",
    reml_bivar_nocove: bool = False,
    reml_bivar_no_constrain: bool = False,
    grm_cutoff: float = -1.0,
    work_dir: str = "",
) -> str:
    """Run bivariate REML analysis to estimate genetic correlation (rg).

    Estimates the genetic correlation between two traits using bivariate
    GREML.  Requires a phenotype file with at least two columns.

    Args:
        grm: GRM file prefix.
        pheno: Phenotype file (with at least 2 trait columns).
        out: Output prefix.  Produces {out}.hsq.
        qcovar: Quantitative covariate file.
        covar: Categorical covariate file.
        mphen: First phenotype column number (default 1).
        mphen2: Second phenotype column number (default 2).
        thread_num: Number of CPU threads.
        reml_maxit: Maximum REML iterations.
        reml_bivar_prevalence: Two prevalence values (space-separated).
        reml_bivar_nocove: If True, ignore residual covariance.
        reml_bivar_no_constrain: If True, no constrain on variance components.
        grm_cutoff: GRM cutoff (-1 = no cutoff).
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--grm", grm, "--pheno", pheno,
            "--thread-num", str(thread_num),
            "--reml-bivar", str(mphen), str(mphen2),
            "--reml-maxit", str(reml_maxit), "--out", out]

    args = _build_args(base, qcovar=qcovar or None, covar=covar or None,
                       reml_bivar_prevalence=reml_bivar_prevalence or None,
                       reml_bivar_nocove=reml_bivar_nocove,
                       reml_bivar_no_constrain=reml_bivar_no_constrain,
                       grm_cutoff=grm_cutoff if grm_cutoff > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_hereg(
    grm: str,
    pheno: str,
    out: str = "gcta_hereg",
    qcovar: str = "",
    covar: str = "",
    mpheno: int = 0,
    thread_num: int = 1,
    grm_cutoff: float = -1.0,
    work_dir: str = "",
) -> str:
    """Run Haseman-Elston regression to estimate variance components.

    HE regression is a method-of-moments estimator that is faster but
    less precise than REML.  Useful as a sanity check or when REML
    fails to converge.

    Args:
        grm: GRM file prefix.
        pheno: Phenotype file path.
        out: Output prefix.
        qcovar: Quantitative covariate file.
        covar: Categorical covariate file.
        mpheno: Phenotype column number (0 = first).
        thread_num: Number of CPU threads.
        grm_cutoff: GRM cutoff (-1 = no cutoff).
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--grm", grm, "--pheno", pheno,
            "--thread-num", str(thread_num), "--HEreg",
            "--out", out]

    args = _build_args(base, qcovar=qcovar or None, covar=covar or None,
                       mpheno=mpheno or None,
                       grm_cutoff=grm_cutoff if grm_cutoff > 0 else None)
    return gcta_run(args, work_dir)


# ============================================================================
# PCA tools
# ============================================================================

@mcp.tool()
def gcta_pca(
    grm: str,
    out: str = "gcta_pca",
    pc_num: int = 20,
    thread_num: int = 1,
    work_dir: str = "",
) -> str:
    """Perform Principal Component Analysis (PCA) from a GRM.

    Computes principal components (PCs) from the genetic relationship
    matrix to capture population structure.  Output eigenvalues and
    eigenvectors are saved in {out}.eigenval and {out}.eigenvec.

    Args:
        grm: GRM file prefix (from gcta_make_grm).
        out: Output file prefix.
        pc_num: Number of principal components to output (default 20).
        thread_num: Number of CPU threads.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    args = (f"--grm {grm} --pca {pc_num} "
            f"--thread-num {thread_num} --out {out}")
    return gcta_run(args, work_dir)


# ============================================================================
# Association analysis tools
# ============================================================================

@mcp.tool()
def gcta_mlma(
    bfile: str,
    grm: str,
    pheno: str,
    out: str = "gcta_mlma",
    qcovar: str = "",
    covar: str = "",
    mpheno: int = 0,
    thread_num: int = 1,
    reml_maxit: int = 100,
    loco: bool = False,
    maf: float = 0.0,
    autosome_num: int = 22,
    work_dir: str = "",
) -> str:
    """Run Mixed Linear Model Association (MLMA) analysis.

    Performs genome-wide association analysis using a linear mixed model
    that accounts for population structure via the GRM.

    Set ``loco=True`` to use the Leave-One-Chromosome-Out (LOCO) method,
    which is more computationally efficient for large datasets.

    Args:
        bfile: PLINK binary file prefix.
        grm: GRM file prefix.
        pheno: Phenotype file path.
        out: Output file prefix.  Produces {out}.mlma (or .loco.mlma).
        qcovar: Quantitative covariate file.
        covar: Categorical covariate file.
        mpheno: Phenotype column number (0 = first).
        thread_num: Number of CPU threads.
        reml_maxit: Maximum REML iterations.
        loco: If True, use MLMA-LOCO method.
        maf: Minor allele frequency filter (0 = no filter).
        autosome_num: Number of autosomes.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    if loco:
        mlma_flag = "--mlma-loco"
    else:
        mlma_flag = "--mlma"

    base = ["--bfile", bfile, "--grm", grm, "--pheno", pheno,
            "--autosome-num", str(autosome_num),
            "--thread-num", str(thread_num),
            mlma_flag, "--reml-maxit", str(reml_maxit), "--out", out]

    args = _build_args(base, qcovar=qcovar or None, covar=covar or None,
                       mpheno=mpheno or None, maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_cojo(
    bfile: str,
    cojo_file: str,
    out: str = "gcta_cojo",
    method: str = "slct",
    cojo_p: float = 5e-8,
    cojo_wind: int = 10000,
    cojo_collinear: float = 0.9,
    cojo_cond_snplist: str = "",
    cojo_gc: float = -1.0,
    cojo_sblup_fac: float = -1.0,
    cojo_top_snps: int = -1,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Run Conditional and Joint (COJO) analysis of GWAS summary data.

    COJO performs stepwise selection, conditional analysis, or joint
    analysis of GWAS summary statistics, using an LD reference panel
    from PLINK bed data.

    Methods:
        "slct"     - Stepwise selection (backward + forward)
        "forward"  - Forward-only selection
        "backward" - Backward-only selection
        "cond"     - Conditional analysis (condition on given SNPs)
        "joint"    - Joint analysis (all SNPs in one model)
        "sblup"    - SBLUP (summary-data BLUP)

    Args:
        bfile: PLINK binary file prefix (LD reference panel).
        cojo_file: GWAS summary data file (.ma format).
        out: Output file prefix.
        method: COJO method (see above).
        cojo_p: P-value threshold for selection (default 5e-8).
        cojo_wind: Window size in Kb for LD calculation (default 10000).
        cojo_collinear: Collinearity threshold (default 0.9).
        cojo_cond_snplist: SNP list file for conditional analysis.
        cojo_gc: Genomic control inflation factor (-1 = not used).
        cojo_sblup_fac: SBLUP factor (for sblup method).
        cojo_top_snps: Number of top SNPs to select (-1 = no limit).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    method = method.lower().strip()
    method_map = {
        "slct": "--cojo-slct",
        "stepwise": "--cojo-slct",
        "forward": "--cojo-forward",
        "backward": "--cojo-backward",
        "joint": "--cojo-joint",
        "sblup": "--cojo-sblup",
    }
    if method not in method_map and method != "cond":
        return json.dumps({
            "error": f"Unknown method '{method}'. Use: slct, forward, backward, cond, joint, sblup"
        }, indent=2)

    base = ["--bfile", bfile, "--cojo-file", cojo_file,
            "--thread-num", str(thread_num), "--out", out]

    if method == "cond":
        base.extend(["--cojo-cond", cojo_cond_snplist])
    elif method == "sblup":
        base.extend([method_map[method], str(cojo_sblup_fac)])
    else:
        base.append(method_map[method])

    args = _build_args(base, cojo_p=cojo_p, cojo_wind=cojo_wind,
                       cojo_collinear=cojo_collinear,
                       cojo_gc=cojo_gc if cojo_gc > 0 else None,
                       cojo_top_snps=cojo_top_snps if cojo_top_snps > 0 else None,
                       maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_gsmr(
    bfile: str,
    gsmr_file: str,
    out: str = "gcta_gsmr",
    gsmr_direction: int = 0,
    gsmr2_beta: bool = False,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Run GSMR (Generalized Summary-data Mendelian Randomization) analysis.

    GSMR estimates the causal effect of an exposure on an outcome using
    GWAS summary data and an LD reference panel.

    Args:
        bfile: PLINK binary file prefix (LD reference).
        gsmr_file: File listing exposure and outcome GWAS summary files.
        out: Output file prefix.
        gsmr_direction: 0 = forward-GSMR, 1 = reverse-GSMR, 2 = bi-GSMR.
        gsmr2_beta: If True, use GSMR2-beta version (multi-SNP HEIDI).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--gsmr-file", gsmr_file,
            "--gsmr-direction", str(gsmr_direction),
            "--thread-num", str(thread_num), "--out", out]

    args = _build_args(base, gsmr2_beta=gsmr2_beta,
                       maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_mtcojo(
    bfile: str,
    mtcojo_file: str,
    mtcojo_bxy: str,
    out: str = "gcta_mtcojo",
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Run mtCOJO (Multi-trait COJO) analysis.

    mtCOJO performs conditional analysis across multiple traits using
    GWAS summary data, adjusting for the genetic correlation between
    traits.

    Args:
        bfile: PLINK binary file prefix (LD reference).
        mtcojo_file: File listing trait GWAS summary files.
        mtcojo_bxy: File with causal effect estimates between traits.
        out: Output file prefix.
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--mtcojo-file", mtcojo_file,
            "--mtcojo-bxy", mtcojo_bxy,
            "--thread-num", str(thread_num), "--out", out]

    args = _build_args(base, maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_fastgwa(
    bfile: str,
    pheno: str,
    out: str = "gcta_fastgwa",
    qcovar: str = "",
    covar: str = "",
    grm: str = "",
    method: str = "fastGWA-mlm",
    thread_num: int = 1,
    maf: float = 0.0,
    autosome_num: int = 22,
    work_dir: str = "",
) -> str:
    """Run fastGWA genome-wide association analysis.

    fastGWA is an efficient mixed-model association method that uses
    sparse GRM techniques for improved speed on large datasets.

    Methods:
        "fastGWA"           - fastGWA (requires --grm)
        "fastGWA-mlm"       - fastGWA-MLM (linear mixed model)
        "fastGWA-mlm-exact" - fastGWA-MLM exact mode
        "fastGWA-lr"        - fastGWA linear regression (no GRM)

    Args:
        bfile: PLINK binary file prefix.
        pheno: Phenotype file path.
        out: Output file prefix.
        qcovar: Quantitative covariate file.
        covar: Categorical covariate file.
        grm: GRM file prefix (for fastGWA with GRM).
        method: fastGWA method (see above).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        autosome_num: Number of autosomes.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    method = method.lower().strip()
    method_map = {
        "fastgwa": "--fastGWA",
        "fastgwa-mlm": "--fastGWA-mlm",
        "fastgwa-mlm-exact": "--fastGWA-mlm-exact",
        "fastgwa-lr": "--fastGWA-lr",
    }
    if method not in method_map:
        return json.dumps({
            "error": f"Unknown method '{method}'. Use: fastGWA, fastGWA-mlm, fastGWA-mlm-exact, fastGWA-lr"
        }, indent=2)

    base = ["--bfile", bfile, "--pheno", pheno,
            "--autosome-num", str(autosome_num),
            "--thread-num", str(thread_num),
            method_map[method], "--out", out]

    args = _build_args(base, qcovar=qcovar or None, covar=covar or None,
                       grm=grm or None, maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_fastbat(
    bfile: str,
    fastbat_file: str,
    out: str = "gcta_fastbat",
    fastbat_gene_list: str = "",
    fastbat_set_list: str = "",
    fastbat_wind: int = 50,
    fastbat_ld_cutoff: float = 0.9,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Run fastBAT gene-based association test.

    fastBAT performs gene-based association tests by combining
    SNP-level p-values within genes, accounting for LD structure.

    Args:
        bfile: PLINK binary file prefix (LD reference).
        fastbat_file: GWAS summary data file (.ma format).
        out: Output file prefix.
        fastbat_gene_list: Gene annotation file (gene -> SNPs mapping).
        fastbat_set_list: SNP set file (predefined sets).
        fastbat_wind: Window size in Kb for gene boundaries (default 50).
        fastbat_ld_cutoff: LD r2 cutoff for removing correlated SNPs (default 0.9).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--fastBAT", fastbat_file,
            "--thread-num", str(thread_num), "--out", out]

    args = _build_args(base, fastBAT_gene_list=fastbat_gene_list or None,
                       fastBAT_set_list=fastbat_set_list or None,
                       fastBAT_wind=fastbat_wind,
                       fastBAT_ld_cutoff=fastbat_ld_cutoff,
                       maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


# ============================================================================
# Simulation tools
# ============================================================================

@mcp.tool()
def gcta_simu_qt(
    bfile: str,
    out: str = "gcta_simu_qt",
    simu_hsq: float = 0.1,
    simu_rep: int = 1,
    simu_causal_loci: str = "",
    simu_seed: float = -1.0,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Simulate quantitative trait phenotypes based on real genotype data.

    Uses the genotype data from PLINK bed files to simulate quantitative
    traits with a specified heritability.  Useful for power analysis and
    method evaluation.

    Args:
        bfile: PLINK binary file prefix.
        out: Output file prefix.
        simu_hsq: Simulated heritability (default 0.1).
        simu_rep: Number of simulation repetitions (default 1).
        simu_causal_loci: File listing causal loci.
        simu_seed: Random seed for simulation (-1 = auto).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--simu-qt",
            "--simu-hsq", str(simu_hsq),
            "--simu-rep", str(simu_rep),
            "--thread-num", str(thread_num), "--out", out]

    args = _build_args(base, simu_causal_loci=simu_causal_loci or None,
                       simu_seed=simu_seed if simu_seed > 0 else None,
                       maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_simu_cc(
    bfile: str,
    simu_case_num: int,
    simu_control_num: int,
    out: str = "gcta_simu_cc",
    simu_hsq: float = 0.1,
    simu_k: float = 0.1,
    simu_rep: int = 1,
    simu_causal_loci: str = "",
    simu_seed: float = -1.0,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Simulate case-control phenotypes based on real genotype data.

    Uses genotype data from PLINK bed files to simulate case-control
    phenotypes with specified heritability and disease prevalence.

    Args:
        bfile: PLINK binary file prefix.
        simu_case_num: Number of cases to simulate.
        simu_control_num: Number of controls to simulate.
        out: Output file prefix.
        simu_hsq: Simulated heritability (default 0.1).
        simu_k: Simulated disease prevalence (default 0.1).
        simu_rep: Number of simulation repetitions (default 1).
        simu_causal_loci: File listing causal loci.
        simu_seed: Random seed (-1 = auto).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--simu-cc",
            str(simu_case_num), str(simu_control_num),
            "--simu-hsq", str(simu_hsq),
            "--simu-k", str(simu_k),
            "--simu-rep", str(simu_rep),
            "--thread-num", str(thread_num), "--out", out]

    args = _build_args(base, simu_causal_loci=simu_causal_loci or None,
                       simu_seed=simu_seed if simu_seed > 0 else None,
                       maf=maf if maf > 0 else None)
    return gcta_run(args, work_dir)


# ============================================================================
# Population genetics tools
# ============================================================================

@mcp.tool()
def gcta_fst(
    bfile: str,
    sub_popu: str,
    out: str = "gcta_fst",
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Calculate Fst (fixation index) for population differentiation analysis.

    Computes pairwise Fst between subpopulations defined in the
    sub-population file.

    Args:
        bfile: PLINK binary file prefix.
        sub_popu: Subpopulation assignment file (FID, IID, pop ID).
        out: Output file prefix.
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    args = (f"--bfile {bfile} --fst --sub-popu {sub_popu} "
            f"--thread-num {thread_num} --out {out}")
    if maf > 0:
        args += f" --maf {maf}"
    return gcta_run(args, work_dir)


# ============================================================================
# Data management tools
# ============================================================================

@mcp.tool()
def gcta_make_bed(
    bfile: str,
    out: str = "gcta_bed",
    maf: float = 0.0,
    thread_num: int = 1,
    extract: str = "",
    exclude: str = "",
    chr_filter: int = 0,
    autosome: bool = False,
    autosome_num: int = 22,
    work_dir: str = "",
) -> str:
    """Convert / filter genotype data and output PLINK bed format.

    Performs data management operations: SNP/individual filtering,
    chromosome extraction, and output in PLINK binary format.

    Args:
        bfile: PLINK binary file prefix (input).
        out: Output file prefix.
        maf: Minor allele frequency filter (0 = no filter).
        thread_num: Number of CPU threads.
        extract: File of SNP IDs to include.
        exclude: File of SNP IDs to exclude.
        chr_filter: Restrict to a specific chromosome (0 = no restriction).
        autosome: If True, restrict to autosomes.
        autosome_num: Number of autosomes.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    base = ["--bfile", bfile, "--make-bed",
            "--autosome-num", str(autosome_num),
            "--thread-num", str(thread_num), "--out", out]

    args = _build_args(base, maf=maf if maf > 0 else None,
                       extract=extract or None, exclude=exclude or None,
                       chr=chr_filter if chr_filter > 0 else None,
                       autosome=autosome)
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_ld_pruning(
    bfile: str,
    ld_pruning_rsq: float = 0.1,
    out: str = "gcta_ld_prune",
    ld_wind: float = 10000.0,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Perform LD pruning on genotype data.

    Removes SNPs in high LD with each other, retaining a set of
    approximately independent SNPs.  Useful for reducing redundancy
    before downstream analyses.

    Args:
        bfile: PLINK binary file prefix.
        ld_pruning_rsq: LD r2 threshold for pruning (default 0.1).
        out: Output file prefix.
        ld_wind: LD window size in Kb (default 10000).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    args = (f"--bfile {bfile} --ld-pruning {ld_pruning_rsq} "
            f"--ld-wind {ld_wind} "
            f"--thread-num {thread_num} --out {out}")
    if maf > 0:
        args += f" --maf {maf}"
    return gcta_run(args, work_dir)


@mcp.tool()
def gcta_ld_score(
    bfile: str,
    out: str = "gcta_ld_score",
    ld_wind: float = 10000.0,
    thread_num: int = 1,
    maf: float = 0.0,
    work_dir: str = "",
) -> str:
    """Calculate LD scores for SNPs in genotype data.

    Computes LD scores (sum of LD r2 with neighboring SNPs) for each
    SNP.  Used as input for LD score regression analyses.

    Args:
        bfile: PLINK binary file prefix.
        out: Output file prefix.
        ld_wind: LD window size in Kb (default 10000).
        thread_num: Number of CPU threads.
        maf: Minor allele frequency filter.
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    args = (f"--bfile {bfile} --ld-score --ld-wind {ld_wind} "
            f"--thread-num {thread_num} --out {out}")
    if maf > 0:
        args += f" --maf {maf}"
    return gcta_run(args, work_dir)


# ============================================================================
# ACAT tools
# ============================================================================

@mcp.tool()
def gcta_acat(
    gene_list: str,
    snp_list: str,
    out: str = "acat_res.csv",
    max_maf: float = 0.01,
    min_mac: int = 20,
    wind: int = 0,
    work_dir: str = "",
) -> str:
    """Run ACAT (Aggregated Cauchy Association Test) gene-based test.

    ACAT combines p-values from multiple SNPs within a gene using a
    Cauchy combination method, providing a gene-level association test.

    Args:
        gene_list: Gene list file (gene -> chromosome, start, end).
        snp_list: SNP list file with p-values.
        out: Output file (default acat_res.csv).
        max_maf: Maximum MAF for included SNPs (default 0.01).
        min_mac: Minimum minor allele count (default 20).
        wind: Extension length for gene boundaries (default 0).
        work_dir: Working directory for GCTA.

    Returns:
        JSON with execution results, log content, and output file list.
    """
    args = (f"--acat --gene-list {gene_list} --snp-list {snp_list} "
            f"--max-maf {max_maf} --min-mac {min_mac} "
            f"--wind {wind} --out {out}")
    return gcta_run(args, work_dir)


# ============================================================================
# File management tools
# ============================================================================

@mcp.tool()
def gcta_list_files(directory: str = ".", pattern: str = "*") -> str:
    """List files in a directory, optionally filtered by pattern.

    Args:
        directory: Directory path (default: current directory).
        pattern: Glob pattern (default: "*" = all files).
                 Examples: "*.log", "*.hsq", "*.grm.*", "test*"

    Returns:
        JSON with list of files and their sizes.
    """
    if not os.path.exists(directory):
        return json.dumps({"error": f"Directory not found: {directory}"}, indent=2)

    files = []
    search_pattern = os.path.join(directory, pattern)
    for f in sorted(glob.glob(search_pattern)):
        if os.path.isfile(f):
            files.append({
                "file": os.path.basename(f),
                "path": os.path.abspath(f),
                "size_bytes": os.path.getsize(f),
            })

    return json.dumps({
        "directory": os.path.abspath(directory),
        "pattern": pattern,
        "file_count": len(files),
        "files": files,
    }, indent=2, ensure_ascii=False)


@mcp.tool()
def gcta_read_file(filepath: str, max_lines: int = 200) -> str:
    """Read and return the content of a text file.

    Useful for inspecting GCTA output files such as .log, .hsq, .mlma,
    .ma, etc.

    Args:
        filepath: Path to the file to read.
        max_lines: Maximum number of lines to return (default 200).

    Returns:
        JSON with file content (truncated to max_lines).
    """
    if not os.path.exists(filepath):
        return json.dumps({"error": f"File not found: {filepath}"}, indent=2)

    try:
        with open(filepath, "r", errors="replace") as f:
            lines = []
            for i, line in enumerate(f):
                if i >= max_lines:
                    lines.append(f"... (truncated at {max_lines} lines)")
                    break
                lines.append(line.rstrip("\n\r"))
        return json.dumps({
            "file": os.path.abspath(filepath),
            "lines_shown": len(lines),
            "content": "\n".join(lines),
        }, indent=2, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e)}, indent=2)


# ============================================================================
# Help documentation
# ============================================================================

def _get_help_docs() -> dict:
    """Return a dictionary of help documentation by topic."""
    return {
        "all": _help_all(),
        "data": _help_data(),
        "grm": _help_grm(),
        "reml": _help_reml(),
        "bivariate_reml": _help_bivariate_reml(),
        "hereg": _help_hereg(),
        "pca": _help_pca(),
        "mlma": _help_mlma(),
        "cojo": _help_cojo(),
        "gsmr": _help_gsmr(),
        "mtcojo": _help_mtcojo(),
        "fastgwa": _help_fastgwa(),
        "fastbat": _help_fastbat(),
        "simulation": _help_simulation(),
        "ld": _help_ld(),
        "fst": _help_fst(),
        "expression": _help_expression(),
        "acat": _help_acat(),
        "options": _help_options(),
    }


def _help_all():
    return """GCTA (Genome-wide Complex Trait Analysis) - Available Analyses

GCTA is a tool for genome-wide association study (GWAS) analyses. It supports:

1. DATA MANAGEMENT (gcta_make_bed, gcta_run)
   - Read PLINK binary PED format (.bed/.bim/.fam)
   - Filter SNPs and individuals
   - Output filtered data in PLINK bed format
   - Calculate allele frequencies, recode genotypes

2. GRM CONSTRUCTION (gcta_make_grm)
   - Calculate Genetic Relationship Matrix (GRM) from genotype data
   - Support for additive, dominance, X-chromosome, and inbreeding GRMs
   - Multi-component GRM analysis

3. REML / GREML (gcta_reml)
   - Estimate SNP-based heritability (h2_SNP) using REML
   - Support for covariates, prevalence, and prior specification
   - BLUP prediction of random effects

4. BIVARIATE REML (gcta_bivariate_reml)
   - Estimate genetic correlation (rg) between two traits
   - Bivariate GREML analysis

5. HE REGRESSION (gcta_hereg)
   - Haseman-Elston regression for variance component estimation
   - Faster but less precise than REML

6. PCA (gcta_pca)
   - Principal component analysis from GRM
   - Population structure analysis

7. MLMA (gcta_mlma)
   - Mixed Linear Model Association analysis
   - MLMA-LOCO (Leave-One-Chromosome-Out) for efficiency
   - Genome-wide association testing

8. COJO (gcta_cojo)
   - Conditional and joint analysis of GWAS summary data
   - Stepwise selection, conditional analysis, joint analysis
   - SBLUP (summary-data BLUP)

9. GSMR (gcta_gsmr)
   - Generalized Summary-data Mendelian Randomization
   - Causal effect estimation using GWAS summary data

10. mtCOJO (gcta_mtcojo)
    - Multi-trait COJO analysis
    - Conditional analysis across multiple traits

11. fastGWA (gcta_fastgwa)
    - Efficient mixed-model association analysis
    - fastGWA-MLM, fastGWA-MLM-exact, fastGWA-LR methods
    - Sparse GRM for improved speed

12. fastBAT (gcta_fastbat)
    - Gene-based association test
    - Combines SNP-level p-values within genes

13. SIMULATION (gcta_simu_qt, gcta_simu_cc)
    - Simulate quantitative trait phenotypes
    - Simulate case-control phenotypes
    - Based on real genotype data

14. LD ANALYSIS (gcta_ld_pruning, gcta_ld_score)
    - LD pruning to identify independent SNPs
    - LD score calculation

15. FST (gcta_fst)
    - Fixation index for population differentiation analysis

16. ACAT (gcta_acat)
    - Aggregated Cauchy Association Test
    - Gene-based test combining SNP p-values

17. EXPRESSION DATA (gcta_run)
    - Expression-based relationship matrix (ERM)
    - ecojo analysis
    - make-erm, make-erm-gz

Use gcta_help(topic) for detailed documentation on any analysis type.
Use gcta_run(args) for full command-line access to all GCTA options.
"""


def _help_data():
    return """GCTA Data Management

INPUT FORMAT:
  --bfile <prefix>     PLINK binary PED format (.bed/.bim/.fam)
  --mbfile <file>      List of multiple bfile prefixes
  --bfile2 <prefix>    Second dataset for comparison
  --dosage-mach <dose> <info>   Mach dosage format
  --dosage-mach-gz <dose> <info>  Mach dosage (gzipped)
  --dosage-beagle <dose> <info>  Beagle dosage format

FILTERING:
  --maf <val>          Minor allele frequency filter (0-0.5)
  --max-maf <val>      Maximum MAF filter
  --chr <num>          Restrict to chromosome
  --autosome           Restrict to autosomes
  --autosome-num <num> Number of autosomes (default 22)
  --extract <file>     Include only SNPs in file
  --exclude <file>     Exclude SNPs in file
  --keep <file>        Keep only individuals in file
  --remove <file>      Remove individuals in file

OUTPUT:
  --make-bed           Output in PLINK bed format
  --freq               Calculate and output allele frequencies
  --recode             Recode genotypes
  --recode-nomiss      Recode without missing data
  --recode-std         Recode with standardization
  --out <prefix>       Output file prefix

EXAMPLES:
  # Filter SNPs by MAF and output PLINK bed
  gcta_run("--bfile test --make-bed --maf 0.05 --out filtered")

  # Extract chromosome 1 only
  gcta_run("--bfile test --make-bed --chr 1 --out chr1")

  # Calculate allele frequencies
  gcta_run("--bfile test --freq --out freq")
"""


def _help_grm():
    return """GCTA GRM (Genetic Relationship Matrix) Construction

The GRM estimates pairwise genomic relatedness between individuals based
on genome-wide SNPs.  It is the foundational step for most downstream
analyses (REML, PCA, MLMA).

BASIC COMMAND:
  --bfile <prefix> --make-grm --out <output_prefix>

OUTPUT FILES:
  {prefix}.grm.bin    - GRM binary data (lower triangle, column-major)
  {prefix}.grm.N.bin  - Number of SNPs used (same format)
  {prefix}.grm.id     - Individual IDs (FID, IID)

OPTIONS:
  --make-grm-alg <0|1>  Algorithm: 0 = VanRaden (default), 1 = alternative
  --grm-adj <val>       GRM adjustment factor (0-1)
  --grm-cutoff <val>    Remove pairs with GRM > cutoff
  --dc <0|1>            Dosage compensation (0=none, 1=X chr)
  --dominance           Calculate dominance GRM
  --make-grm-xchr       Calculate X-chromosome GRM
  --make-grm-inbred     Calculate inbreeding GRM
  --save-ram            Use memory-saving mode
  --autosome-num <num>  Number of autosomes
  --autosome            Restrict to autosomes
  --chr <num>           Restrict to specific chromosome
  --maf <val>           MAF filter
  --thread-num <num>    Number of CPU threads

MULTI-GRM:
  --mgrm <file>         Multiple GRM file list
  --mgrm-bin <file>     Multiple GRM binary file list
  --mgrm-gz <file>      Multiple GRM gzipped file list

EXAMPLES:
  # Basic GRM
  gcta_make_grm(bfile="test", out="test_grm")

  # GRM with MAF filter and 4 threads
  gcta_make_grm(bfile="test", out="test_grm", maf=0.01, thread_num=4)

  # Dominance GRM
  gcta_make_grm(bfile="test", out="test_dgrm", dominance=True)
"""


def _help_reml():
    return """GCTA REML / GREML Analysis

REML (Restricted Maximum Likelihood) estimates the proportion of
phenotypic variance explained by all genome-wide SNPs (SNP-based
heritability, h2_SNP).

BASIC COMMAND:
  --grm <prefix> --pheno <file> --reml --out <output_prefix>

REQUIRED INPUT:
  --grm <prefix>        GRM file prefix (from gcta_make_grm)
  --pheno <file>        Phenotype file (FID, IID, phenotype)
  --reml                Perform REML analysis
  --out <prefix>        Output file prefix

OUTPUT FILES:
  {prefix}.hsq         - Heritability estimates and SE
  {prefix}.log         - Analysis log

COVARIATES:
  --qcovar <file>      Quantitative covariate file
  --covar <file>       Categorical covariate file

REML OPTIONS:
  --reml-maxit <num>   Maximum iterations (default 100)
  --reml-priors <vals> Prior values of variance explained
  --reml-priors-var <vals>  Prior values of variance components
  --reml-no-constrain  Do not constrain variance components
  --reml-no-lrt        Skip likelihood ratio test
  --reml-pred-rand     Predict random effects (BLUP)
  --reml-est-fix       Estimate fixed effects

OTHER OPTIONS:
  --prevalence <val>   Disease prevalence for liability conversion
  --grm-cutoff <val>   GRM cutoff for related pairs
  --gxe <file>         GxE interaction file
  --mpheno <num>       Phenotype column number
  --thread-num <num>   Number of CPU threads

EXAMPLES:
  # Basic REML
  gcta_reml(grm="test_grm", pheno="test.phen", out="test_reml")

  # REML with covariates
  gcta_reml(grm="test_grm", pheno="test.phen", out="test_reml",
            qcovar="covariates.txt", covar="sex.txt")
"""


def _help_bivariate_reml():
    return """GCTA Bivariate REML Analysis

Bivariate REML estimates the genetic correlation (rg) between two traits
using bivariate GREML.

BASIC COMMAND:
  --grm <prefix> --pheno <file> --reml-bivar <mphen> <mphen2> --out <prefix>

REQUIRED INPUT:
  --grm <prefix>            GRM file prefix
  --pheno <file>            Phenotype file (with 2+ trait columns)
  --reml-bivar <p1> <p2>    Bivariate REML with trait columns p1 and p2
  --out <prefix>            Output file prefix

OUTPUT FILES:
  {prefix}.hsq             - Variance components and genetic correlation

OPTIONS:
  --reml-maxit <num>       Maximum iterations
  --reml-bivar-prevalence <k1> <k2>  Prevalence for both traits
  --reml-bivar-nocove      Ignore residual covariance
  --reml-bivar-no-constrain  No constrain on variance components
  --qcovar <file>          Quantitative covariate file
  --covar <file>           Categorical covariate file
  --grm-cutoff <val>       GRM cutoff
  --thread-num <num>       Number of CPU threads

EXAMPLE:
  gcta_bivariate_reml(grm="test_grm", pheno="traits.phen",
                      out="bireml", mphen=1, mphen2=2)
"""


def _help_hereg():
    return """GCTA Haseman-Elston (HE) Regression

HE regression is a method-of-moments estimator for variance components.
It is faster but less precise than REML.  Useful as a sanity check or
when REML fails to converge.

BASIC COMMAND:
  --grm <prefix> --pheno <file> --HEreg --out <prefix>

OUTPUT FILES:
  {prefix}.hsq    - HE regression estimates

OPTIONS:
  --qcovar <file>        Quantitative covariate file
  --covar <file>         Categorical covariate file
  --mpheno <num>         Phenotype column number
  --grm-cutoff <val>     GRM cutoff
  --thread-num <num>     Number of CPU threads

EXAMPLE:
  gcta_hereg(grm="test_grm", pheno="test.phen", out="test_he")
"""


def _help_pca():
    return """GCTA Principal Component Analysis (PCA)

PCA from the GRM captures population structure.  Output eigenvalues
and eigenvectors are saved.

BASIC COMMAND:
  --grm <prefix> --pca <num_pcs> --out <prefix>

OUTPUT FILES:
  {prefix}.eigenval    - Eigenvalues (one per line)
  {prefix}.eigenvec    - Eigenvectors (FID, IID, PC1, PC2, ...)

OPTIONS:
  --pca <num>          Number of PCs (default 20)
  --thread-num <num>   Number of CPU threads

EXAMPLE:
  gcta_pca(grm="test_grm", out="test_pca", pc_num=10)
"""


def _help_mlma():
    return """GCTA Mixed Linear Model Association (MLMA)

MLMA performs genome-wide association analysis using a linear mixed
model that accounts for population structure via the GRM.

METHODS:
  --mlma           Standard MLMA (uses full GRM)
  --mlma-loco      MLMA with Leave-One-Chromosome-Out (faster, recommended)

BASIC COMMAND:
  --bfile <prefix> --grm <prefix> --pheno <file> --mlma-loco --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>     PLINK binary file prefix (genotype data)
  --grm <prefix>       GRM file prefix
  --pheno <file>       Phenotype file
  --out <prefix>       Output file prefix

OUTPUT FILES:
  {prefix}.mlma         - Association results (Chr, SNP, BP, A1, F, BETA, SE, P)

COVARIATES:
  --qcovar <file>      Quantitative covariate file
  --covar <file>       Categorical covariate file

OPTIONS:
  --mpheno <num>       Phenotype column number
  --reml-maxit <num>   Maximum REML iterations
  --maf <val>          Minor allele frequency filter
  --autosome-num <num> Number of autosomes
  --thread-num <num>   Number of CPU threads

EXAMPLES:
  # MLMA-LOCO (recommended)
  gcta_mlma(bfile="test", grm="test_grm", pheno="test.phen",
            out="test_mlma", loco=True)

  # Standard MLMA
  gcta_mlma(bfile="test", grm="test_grm", pheno="test.phen",
            out="test_mlma", loco=False)
"""


def _help_cojo():
    return """GCTA Conditional and Joint (COJO) Analysis

COJO performs stepwise selection, conditional analysis, or joint analysis
of GWAS summary statistics, using an LD reference panel from PLINK bed
data.

BASIC COMMAND:
  --bfile <prefix> --cojo-file <file> --cojo-slct --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>         PLINK binary file prefix (LD reference panel)
  --cojo-file <file>       GWAS summary data file (.ma format)
  --out <prefix>           Output file prefix

GWAS SUMMARY FILE FORMAT (.ma):
  Columns: SNP A1 A2 freq b se p n
  (Header line required)

METHODS:
  --cojo-slct              Stepwise selection (backward + forward)
  --cojo-forward           Forward-only selection
  --cojo-backward          Backward-only selection
  --cojo-cond <file>       Conditional analysis (condition on given SNPs)
  --cojo-joint             Joint analysis (all SNPs in one model)
  --cojo-sblup <factor>    SBLUP (summary-data BLUP)

OPTIONS:
  --cojo-p <val>           P-value threshold (default 5e-8)
  --cojo-wind <kb>         Window size in Kb (default 10000)
  --cojo-collinear <val>   Collinearity threshold (default 0.9)
  --cojo-gc <val>          Genomic control inflation factor
  --cojo-top-SNPs <num>    Number of top SNPs to select
  --maf <val>              Minor allele frequency filter
  --thread-num <num>       Number of CPU threads

OUTPUT FILES:
  {prefix}.jma.cojo    - Joint analysis results
  {prefix}.cma.cojo    - Conditional analysis results
  {prefix}.log         - Analysis log

EXAMPLES:
  # Stepwise selection
  gcta_cojo(bfile="test", cojo_file="gwas.ma", out="cojo_res",
            method="slct")

  # Conditional analysis
  gcta_cojo(bfile="test", cojo_file="gwas.ma", out="cojo_res",
            method="cond", cojo_cond_snplist="cond_snps.txt")
"""


def _help_gsmr():
    return """GCTA GSMR (Generalized Summary-data Mendelian Randomization)

GSMR estimates the causal effect of an exposure on an outcome using
GWAS summary data and an LD reference panel.

BASIC COMMAND:
  --bfile <prefix> --gsmr-file <file> --gsmr-direction <0|1|2> --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>             PLINK binary file prefix (LD reference)
  --gsmr-file <file>           File listing exposure and outcome GWAS files
  --gsmr-direction <num>       0 = forward, 1 = reverse, 2 = bi-directional
  --out <prefix>               Output file prefix

GSMR FILE FORMAT:
  Exposure_GWAS_filename  Outcome_GWAS_filename  n_exp  n_out

OPTIONS:
  --gsmr2-beta                 Use GSMR2-beta version (multi-SNP HEIDI)
  --gsmr-snp-min <num>         Minimum number of SNP instruments
  --gwas-thresh <val>          GWAS p-value threshold
  --clump-kb <val>             Clumping window size
  --clump-r2 <val>             Clumping r2 threshold
  --heidi-thresh <val>         HEIDI-outlier p-value threshold
  --diff-freq <val>            Frequency difference threshold
  --thread-num <num>           Number of CPU threads

OUTPUT FILES:
  {prefix}.gsmr         - GSMR results
  {prefix}.log          - Analysis log

EXAMPLE:
  gcta_gsmr(bfile="test", gsmr_file="gsmr_data.txt",
            out="gsmr_res", gsmr_direction=0)
"""


def _help_mtcojo():
    return """GCTA mtCOJO (Multi-trait COJO) Analysis

mtCOJO performs conditional analysis across multiple traits using GWAS
summary data, adjusting for the genetic correlation between traits.

BASIC COMMAND:
  --bfile <prefix> --mtcojo-file <file> --mtcojo-bxy <file> --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>         PLINK binary file prefix (LD reference)
  --mtcojo-file <file>     File listing trait GWAS summary files
  --mtcojo-bxy <file>      File with causal effect estimates between traits
  --out <prefix>           Output file prefix

OPTIONS:
  --gwas-thresh <val>      GWAS p-value threshold
  --clump-kb <val>         Clumping window size
  --clump-r2 <val>         Clumping r2 threshold
  --diff-freq <val>        Frequency difference threshold
  --heidi-thresh <val>     HEIDI-outlier threshold
  --thread-num <num>       Number of CPU threads

EXAMPLE:
  gcta_mtcojo(bfile="test", mtcojo_file="traits.txt",
              mtcojo_bxy="bxy.txt", out="mtcojo_res")
"""


def _help_fastgwa():
    return """GCTA fastGWA Analysis

fastGWA is an efficient mixed-model association method that uses sparse
GRM techniques for improved speed on large datasets.

METHODS:
  --fastGWA           fastGWA (requires --grm sparse GRM)
  --fastGWA-mlm       fastGWA-MLM (linear mixed model, recommended)
  --fastGWA-mlm-exact fastGWA-MLM exact mode
  --fastGWA-lr        fastGWA linear regression (no GRM)

BASIC COMMAND:
  --bfile <prefix> --pheno <file> --fastGWA-mlm --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>     PLINK binary file prefix
  --pheno <file>       Phenotype file
  --out <prefix>       Output file prefix

COVARIATES:
  --qcovar <file>      Quantitative covariate file
  --covar <file>       Categorical covariate file

OPTIONS:
  --grm <prefix>       GRM file prefix (for fastGWA with GRM)
  --maf <val>          Minor allele frequency filter
  --autosome-num <num> Number of autosomes
  --thread-num <num>   Number of CPU threads

OUTPUT FILES:
  {prefix}.fastGWA     - Association results
  {prefix}.log         - Analysis log

EXAMPLES:
  # fastGWA-MLM (recommended)
  gcta_fastgwa(bfile="test", pheno="test.phen", out="fastgwa_res",
               method="fastGWA-mlm")

  # fastGWA-LR (no GRM, fastest)
  gcta_fastgwa(bfile="test", pheno="test.phen", out="fastgwa_res",
               method="fastGWA-lr")
"""


def _help_fastbat():
    return """GCTA fastBAT Gene-based Association Test

fastBAT performs gene-based association tests by combining SNP-level
p-values within genes, accounting for LD structure.

BASIC COMMAND:
  --bfile <prefix> --fastBAT <gwas_summary> --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>               PLINK binary file prefix (LD reference)
  --fastBAT <file>               GWAS summary data file (.ma format)
  --fastBAT-gene-list <file>     Gene annotation file
  --fastBAT-set-list <file>      SNP set file (predefined sets)
  --out <prefix>                 Output file prefix

GWAS SUMMARY FILE FORMAT (.ma):
  Columns: SNP A1 A2 freq b se p n

OPTIONS:
  --fastBAT-wind <kb>            Window size in Kb (default 50)
  --fastBAT-ld-cutoff <val>      LD r2 cutoff (default 0.9)
  --fastBAT-write-snpset         Write SNP set file
  --maf <val>                    Minor allele frequency filter
  --thread-num <num>             Number of CPU threads

OUTPUT FILES:
  {prefix}.fastBAT    - Gene-based test results
  {prefix}.log        - Analysis log

EXAMPLE:
  gcta_fastbat(bfile="test", fastbat_file="gwas.ma",
               out="fastbat_res", fastbat_gene_list="genes.txt")
"""


def _help_simulation():
    return """GCTA Phenotype Simulation

GCTA can simulate phenotypes based on real genotype data, useful for
power analysis and method evaluation.

QUANTITATIVE TRAIT SIMULATION:
  --bfile <prefix> --simu-qt --simu-hsq <h2> --out <prefix>

CASE-CONTROL SIMULATION:
  --bfile <prefix> --simu-cc <cases> <controls> --simu-hsq <h2>
            --simu-k <prevalence> --out <prefix>

OPTIONS:
  --simu-hsq <val>         Simulated heritability (default 0.1)
  --simu-k <val>           Disease prevalence (default 0.1)
  --simu-rep <num>         Number of repetitions (default 1)
  --simu-causal-loci <file>  File of causal loci
  --simu-seed <val>        Random seed
  --simu-eff-mod <0|1>     Effect model: 0 = additive, 1 = non-additive
  --thread-num <num>       Number of CPU threads

EXAMPLES:
  # Simulate QT with h2=0.5
  gcta_simu_qt(bfile="test", out="simu_qt", simu_hsq=0.5)

  # Simulate CC with 1000 cases, 1000 controls
  gcta_simu_cc(bfile="test", simu_case_num=1000, simu_control_num=1000,
               out="simu_cc", simu_hsq=0.3, simu_k=0.05)
"""


def _help_ld():
    return """GCTA LD Analysis

GCTA provides several LD analysis tools:

LD PRUNING:
  --bfile <prefix> --ld-pruning <r2> --ld-wind <kb> --out <prefix>
  Removes SNPs in high LD, retaining approximately independent SNPs.

LD SCORE CALCULATION:
  --bfile <prefix> --ld-score --ld-wind <kb> --out <prefix>
  Computes LD scores (sum of LD r2 with neighboring SNPs).

LD SCORE REGRESSION:
  --bfile <prefix> --ld-score-region --ld-wind <kb> --out <prefix>
  Performs LD score regression to estimate heritability and confounding.

OPTIONS:
  --ld-pruning <r2>      LD r2 threshold for pruning
  --ld-wind <kb>         LD window size in Kb
  --ld-score             Calculate LD scores
  --ld-score-region      LD score regression
  --ld-score-multi <f>   Multi-set LD score
  --ld-rsq-cutoff <val>  LD r2 cutoff
  --ld-max-rsq           Maximum LD r2
  --ld-step <num>        LD step size
  --thread-num <num>     Number of CPU threads

EXAMPLES:
  # LD pruning with r2=0.1
  gcta_ld_pruning(bfile="test", ld_pruning_rsq=0.1, out="pruned")

  # LD score calculation
  gcta_ld_score(bfile="test", out="ldscores")
"""


def _help_fst():
    return """GCTA Fst (Fixation Index) Analysis

Fst measures population differentiation.  GCTA computes pairwise Fst
between subpopulations.

BASIC COMMAND:
  --bfile <prefix> --fst --sub-popu <file> --out <prefix>

REQUIRED INPUT:
  --bfile <prefix>       PLINK binary file prefix
  --sub-popu <file>      Subpopulation file (FID, IID, pop ID)
  --out <prefix>         Output file prefix

OPTIONS:
  --maf <val>            Minor allele frequency filter
  --thread-num <num>     Number of CPU threads

OUTPUT FILES:
  {prefix}.Fst    - Fst results
  {prefix}.log    - Analysis log

EXAMPLE:
  gcta_fst(bfile="test", sub_popu="pops.txt", out="fst_res")
"""


def _help_expression():
    return """GCTA Expression Data Analysis

GCTA supports analysis of gene expression data:

EXPRESSION FILE INPUT:
  --efile <file>         Expression data file

EXPRESSION RELATIONSHIP MATRIX (ERM):
  --make-erm             Make ERM (binary output)
  --make-erm-gz          Make ERM (gzipped text output)
  --make-erm-alg <1|2|3> Algorithm for ERM calculation

ECOJO ANALYSIS:
  --ecojo <file>         ecojo analysis with MA file
  --ecojo-slct           ecojo stepwise selection
  --ecojo-p <val>        ecojo p-value threshold
  --ecojo-collinear <v>  ecojo collinearity threshold
  --ecojo-blup <lambda>  ecojo BLUP with lambda

E-COR (Expression Correlation):
  --e-cor <file>         Expression correlation file

EXAMPLE:
  gcta_run("--efile expr.txt --make-erm --out expr_erm")
"""


def _help_acat():
    return """GCTA ACAT (Aggregated Cauchy Association Test)

ACAT combines p-values from multiple SNPs within a gene using a Cauchy
combination method, providing a gene-level association test.

BASIC COMMAND:
  --acat --gene-list <file> --snp-list <file> --out <file>

REQUIRED INPUT:
  --gene-list <file>     Gene list file (gene, chr, start, end)
  --snp-list <file>      SNP list file with p-values
  --out <file>           Output file (default acat_res.csv)

OPTIONS:
  --max-maf <val>        Maximum MAF for included SNPs (default 0.01)
  --min-mac <num>        Minimum minor allele count (default 20)
  --wind <num>           Extension length for gene boundaries

OUTPUT FILE:
  {out}                  - ACAT gene-based test results

EXAMPLE:
  gcta_acat(gene_list="genes.txt", snp_list="snps.txt",
            out="acat_results.csv")
"""


def _help_options():
    return """GCTA Complete Command-Line Options Reference

DATA INPUT:
  --bfile <prefix>           PLINK binary PED format
  --mbfile <file>            Multiple bfile list
  --bfile2 <prefix>          Second bfile
  --dosage-mach <d> <i>      Mach dosage
  --dosage-mach-gz <d> <i>   Mach dosage (gzipped)
  --dosage-beagle <d> <i>    Beagle dosage

DATA MANAGEMENT:
  --make-bed                 Output PLINK bed format
  --freq / --freqx           Calculate allele frequencies
  --recode / --recode-nomiss / --recode-std  Recode genotypes
  --keep <file>              Keep individuals
  --remove <file>            Remove individuals
  --extract <file>           Include SNPs
  --exclude <file>           Exclude SNPs
  --chr <num>                Chromosome filter
  --autosome                 Autosome filter
  --autosome-num <num>       Number of autosomes
  --maf <val>                MAF filter
  --max-maf <val>            Max MAF filter
  --update-sex <file>        Update sex
  --update-freq <file>       Update allele frequencies
  --update-ref-allele <f>    Update reference allele
  --save-ram                 Memory-saving mode

GRM:
  --make-grm                 Make GRM (binary output)
  --make-grm-gz              Make GRM (text output)
  --make-grm-alg <0|1>       GRM algorithm
  --make-grm-d               Dominance GRM
  --make-grm-xchr            X-chromosome GRM
  --make-grm-inbred          Inbreeding GRM
  --grm <prefix>             Read GRM
  --grm-gz <prefix>          Read GRM (gzipped)
  --mgrm <file>              Multiple GRMs
  --mgrm-gz <file>           Multiple GRMs (gzipped)
  --grm-cutoff <val>         GRM cutoff
  --grm-adj <val>            GRM adjustment
  --dc <0|1>                 Dosage compensation

PCA:
  --pca <num>                PCA with num PCs
  --pc-loading <file>        PC loading
  --project-loading <f> <n>  Project loading

REML:
  --reml                     REML analysis
  --reml-bivar <p1> <p2>     Bivariate REML
  --reml-maxit <num>         Max iterations
  --reml-priors <vals>       Prior variance explained
  --reml-priors-var <vals>   Prior variance components
  --reml-no-constrain        No constrain
  --reml-no-lrt              No LRT
  --reml-pred-rand           Predict random effects
  --reml-est-fix             Estimate fixed effects
  --prevalence <val>         Disease prevalence

HE REGRESSION:
  --HEreg                    HE regression
  --HEreg-bivar <p1> <p2>    Bivariate HE regression

MLMA:
  --mlma                     MLMA association
  --mlma-loco                MLMA-LOCO
  --mlma-no-adj-covar        No adj covar

COJO:
  --cojo-file <file>         GWAS summary data
  --cojo-slct                Stepwise selection
  --cojo-forward             Forward selection
  --cojo-backward            Backward selection
  --cojo-joint               Joint analysis
  --cojo-cond <file>         Conditional analysis
  --cojo-sblup <factor>      SBLUP
  --cojo-p <val>             P-value threshold
  --cojo-wind <kb>           Window size
  --cojo-collinear <val>     Collinearity threshold
  --cojo-gc <val>            Genomic control
  --cojo-top-SNPs <num>      Top SNPs

GSMR:
  --gsmr-file <file>         GSMR data file
  --gsmr-direction <0|1|2>   Direction
  --gsmr2-beta               GSMR2-beta
  --gsmr-snp-min <num>       Min SNP instruments
  --gwas-thresh <val>        GWAS threshold
  --heidi-thresh <val>       HEIDI threshold

mtCOJO:
  --mtcojo-file <file>       mtCOJO data file
  --mtcojo-bxy <file>        Bxy file

fastGWA:
  --fastGWA                  fastGWA
  --fastGWA-mlm              fastGWA-MLM
  --fastGWA-mlm-exact        fastGWA-MLM exact
  --fastGWA-lr               fastGWA-LR

fastBAT:
  --fastBAT <file>           fastBAT GWAS summary
  --fastBAT-gene-list <f>    Gene list
  --fastBAT-set-list <f>     SNP set list
  --fastBAT-wind <kb>        Window size
  --fastBAT-ld-cutoff <v>    LD cutoff

SIMULATION:
  --simu-qt                   Simulate QT
  --simu-cc <cases> <ctrl>    Simulate CC
  --simu-hsq <val>            Heritability
  --simu-k <val>              Prevalence
  --simu-rep <num>            Repetitions
  --simu-causal-loci <file>   Causal loci
  --simu-seed <val>           Seed

LD:
  --ld <file>                 LD analysis
  --ld-pruning <r2>           LD pruning
  --ld-score                  LD score
  --ld-score-region           LD score regression
  --ld-wind <kb>              LD window
  --ld-rsq-cutoff <val>       LD r2 cutoff

FST:
  --fst                       Fst analysis
  --sub-popu <file>           Subpopulation file

EXPRESSION:
  --efile <file>              Expression file
  --e-cor <file>              Expression correlation
  --make-erm                  Make ERM
  --ecojo <file>              ecojo analysis

ACAT:
  --acat                      ACAT analysis
  --gene-list <file>          Gene list
  --snp-list <file>           SNP list
  --max-maf <val>             Max MAF
  --min-mac <num>             Min MAC

THREADING:
  --thread-num <num>          Number of threads
  --threads <num>             Alias for --thread-num

OUTPUT:
  --out <prefix>              Output prefix
"""


# ============================================================================
# Main entry point
# ============================================================================

def main():
    print("GCTA MCP Server starting...", file=sys.stderr)
    print(f"  Binary path:    {GCTA_BINARY_PATH}", file=sys.stderr)
    print(f"  Work directory: {GCTA_WORK_DIR}", file=sys.stderr)
    print(f"  Timeout:        {GCTA_TIMEOUT}s", file=sys.stderr)
    mcp.run()


if __name__ == "__main__":
    main()
