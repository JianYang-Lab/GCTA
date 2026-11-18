# GCTA MCP Server

Expose all GCTA functionality to AI tools (Claude, Codex, opencode, etc.) via the MCP (Model Context Protocol) protocol, enabling users to invoke various GCTA analysis functions through natural language conversations.

## Directory Structure and Defaults

In the source tree, the server is under `.opencode/gcta-mcp/` and uses the
CMake Release binary and source root by default:

```
GCTA/
├── build/Release/gcta64
├── .opencode/gcta-mcp/
│   ├── server.py
│   ├── pyproject.toml
│   ├── requirements.txt
│   └── README.md
└── tests/data/gcta_legacy/
```

A release package can instead keep `gcta` and `gcta-mcp/` in the same
directory. In that layout, the defaults remain `../gcta` relative to
`server.py` and the package root as the working directory.

`GCTA_BINARY_PATH` and `GCTA_WORK_DIR` override both layouts. The server does
not build GCTA; build `gcta64` first or point it at an existing binary.

## Installation

### 1. Install Python dependencies

```bash
python3 -m pip install -r .opencode/gcta-mcp/requirements.txt
```

### 2. Build GCTA and make the binary executable (if needed)

Follow the root README's CMake Release build instructions, then run:

```bash
chmod +x build/Release/gcta64
```

### 3. Verify

```bash
python3 .opencode/gcta-mcp/server.py
```

If you see the "GCTA MCP Server starting..." message, the server has started successfully.

## Configuration in AI Tools

Replace `<GCTA_DIR>` below with the full path to the `GCTA` folder.

### Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or the corresponding config file for your platform:

```json
{
  "mcpServers": {
    "gcta": {
      "command": "python3",
      "args": ["<GCTA_DIR>/.opencode/gcta-mcp/server.py"]
    }
  }
}
```

### opencode

Edit `opencode.json`:

```json
{
  "mcp": {
    "gcta": {
      "type": "local",
      "command": ["python3", "<GCTA_DIR>/.opencode/gcta-mcp/server.py"],
      "enabled": true
    }
  }
}
```

### Codex (OpenAI)

Edit `~/.codex/config.json`:

```json
{
  "mcpServers": {
    "gcta": {
      "command": "python3",
      "args": ["<GCTA_DIR>/.opencode/gcta-mcp/server.py"]
    }
  }
}
```

## Environment Variables (Optional)

| Variable | Description | Default |
|----------|-------------|---------|
| `GCTA_BINARY_PATH` | Path to the GCTA binary | `build/Release/gcta64` in the source tree; `../gcta` in a release package |
| `GCTA_WORK_DIR` | Working directory (where input/output files reside) | Source or release package root |
| `GCTA_TIMEOUT` | Execution timeout (seconds) | `3600` |

## Tool List

This MCP service provides 23 tools covering all major GCTA analysis functions:

| Tool | Function |
|------|----------|
| `gcta_run` | Generic execution (arbitrary CLI arguments) |
| `gcta_help` | Get GCTA documentation |
| `gcta_check` | Check GCTA configuration |
| `gcta_make_grm` | Build GRM (Genetic Relationship Matrix) |
| `gcta_reml` | REML analysis (estimate variance components/heritability) |
| `gcta_bivariate_reml` | Bivariate REML (estimate genetic correlation) |
| `gcta_hereg` | HE regression analysis |
| `gcta_pca` | Principal component analysis |
| `gcta_mlma` | Mixed linear model association analysis |
| `gcta_cojo` | Conditional and joint analysis |
| `gcta_gsmr` | GSMR Mendelian randomization |
| `gcta_mtcojo` | Multi-trait COJO analysis |
| `gcta_fastgwa` | fastGWA genome-wide association analysis |
| `gcta_fastbat` | fastBAT gene-level association test |
| `gcta_simu_qt` | Simulate quantitative traits |
| `gcta_simu_cc` | Simulate case-control phenotypes |
| `gcta_fst` | Fst population differentiation analysis |
| `gcta_make_bed` | Data management (format conversion, filtering) |
| `gcta_ld_pruning` | LD pruning |
| `gcta_ld_score` | LD score calculation |
| `gcta_acat` | ACAT gene-level test |
| `gcta_list_files` | List files in working directory |
| `gcta_read_file` | Read file contents |

## Usage Examples

Once configured, you can use natural language requests in AI tools:

> "Build a GRM using the test dataset with output prefix test_grm"

From the source root, the AI can call
`gcta_make_grm(bfile="tests/data/gcta_legacy/test", out="test_grm")`.

> "Perform REML analysis on test_grm with phenotype file tests/data/gcta_legacy/test.phen"

The AI will call
`gcta_reml(grm="test_grm", pheno="tests/data/gcta_legacy/test.phen", out="test_reml")`.

> "Read the contents of test_reml.hsq"

The AI will call `gcta_read_file(filepath="test_reml.hsq")`.

## Technical Details

- The source-tree default is the CMake target `gcta64`; its linkage and runtime
  requirements are determined by the selected GCTA build configuration and
  dependencies, not by this MCP service.
- The MCP service uses the Python `mcp` 1.x FastMCP API; both dependency
  files constrain the package to `mcp>=1.0.0,<2.0.0`.
- Relative input and output paths are resolved from `GCTA_WORK_DIR`.
- `pyproject.toml` exposes the same server through the `gcta-mcp` console
  script after installation; direct `python3 server.py` execution is also
  supported. Automatic source/release path detection applies to those two
  on-disk layouts, so set `GCTA_BINARY_PATH` and, when needed,
  `GCTA_WORK_DIR` for an installed console script.

## License

Refer to the GCTA root README and the license notices in the corresponding
source files. This MCP documentation does not alter those licenses.
