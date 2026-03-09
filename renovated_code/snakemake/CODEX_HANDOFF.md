# xQTL Pipeline — Codex Handoff

## Context

This repository (`xqtl-protocol`) implements a bioinformatics pipeline for quantitative trait loci (xQTL) analysis. It has two parallel representations:

1. **Mini-protocols** — top-level Jupyter notebooks in `code/` that document each analysis module step-by-step. Each notebook has numbered sections (i, ii, iii…), and each section calls one `sos run` command for one SoS workflow step. These are the **ground truth** for what steps exist and how they should be parameterized.

2. **Snakemake pipeline** in `renovated_code/snakemake/` — the production pipeline that orchestrates execution. Rules call SoS notebooks via `sos run route3/notebooks/<notebook>.ipynb <step>`.

The `route3/` sub-pipeline is a refactored version where SoS task blocks call external modular scripts (in `renovated_code/`) instead of inline R/Python. `transform_notebooks.py` generates the route3 notebooks from the originals.

---

## Directory Map

```
xqtl-protocol/
├── code/                                      # Mini-protocol notebooks (ground truth)
│   ├── molecular_phenotypes/
│   │   ├── bulk_expression.ipynb              # Mini-protocol: RNA-seq expression
│   │   ├── calling/RNA_calling.ipynb          # Implementation notebook
│   │   └── QC/bulk_expression_QC.ipynb        # Implementation notebook
│   ├── data_preprocessing/
│   │   ├── genotype_preprocessing.ipynb       # Mini-protocol: genotype preprocessing
│   │   ├── phenotype_preprocessing.ipynb      # Mini-protocol: phenotype preprocessing
│   │   ├── covariate_preprocessing.ipynb      # Mini-protocol: covariate preprocessing
│   │   ├── genotype/{VCF_QC,PCA,GWAS_QC,genotype_formatting}.ipynb
│   │   ├── covariate/{covariate_formatting,covariate_hidden_factor}.ipynb
│   │   └── phenotype/{phenotype_formatting,phenotype_imputation}.ipynb
│   ├── association_scan/
│   │   └── qtl_association_testing.ipynb      # Mini-protocol: association testing
│   └── mnm_analysis/
│       └── mnm_miniprotocol.ipynb             # Mini-protocol: fine-mapping
├── pipeline/                                   # Symlinks to code/ notebooks (sos run targets)
├── renovated_code/
│   ├── snakemake/
│   │   ├── Snakefile                           # Main pipeline entry
│   │   ├── config.yaml                         # Pipeline config
│   │   ├── rules/sos/                          # Snakemake rule modules
│   │   │   ├── 01_molecular_phenotypes.smk     # ↔ bulk_expression.ipynb
│   │   │   ├── 02_genotype_preprocessing.smk   # ↔ genotype_preprocessing.ipynb
│   │   │   ├── 03_sample_qc_pca.smk            # ↔ genotype_preprocessing.ipynb (PCA section)
│   │   │   ├── 04_phenotype_covariate_prep.smk # ↔ phenotype_preprocessing.ipynb + covariate_preprocessing.ipynb
│   │   │   ├── 05_association_testing.smk      # ↔ qtl_association_testing.ipynb
│   │   │   └── 06_univariate_finemapping.smk   # ↔ mnm_miniprotocol.ipynb
│   │   ├── route3/
│   │   │   ├── Snakefile                       # Route3 entry point
│   │   │   ├── transform_notebooks.py          # Generates route3/notebooks/ from pipeline/
│   │   │   └── notebooks/                      # Generated wrapper notebooks (sos run these)
│   │   │       ├── RNA_calling.ipynb
│   │   │       ├── bulk_expression_QC.ipynb
│   │   │       ├── bulk_expression_normalization.ipynb
│   │   │       ├── VCF_QC.ipynb
│   │   │       ├── GWAS_QC.ipynb
│   │   │       ├── genotype_formatting.ipynb
│   │   │       ├── PCA.ipynb
│   │   │       ├── phenotype_formatting.ipynb
│   │   │       ├── phenotype_imputation.ipynb
│   │   │       ├── covariate_formatting.ipynb
│   │   │       ├── covariate_hidden_factor.ipynb
│   │   │       ├── TensorQTL.ipynb
│   │   │       ├── mnm_regression.ipynb
│   │   │       └── rss_analysis.ipynb
│   │   └── dryrun/
│   │       ├── mwe/     # 21 dry-run test scripts (NN_<step>.sh)
│   │       └── output/  # Expected output for each test (NN_<step>.txt)
│   └── <category>/      # Modular R/Python/bash scripts called by route3 notebooks
│       ├── calling/RNA_calling.sh
│       ├── genotype/PCA.R
│       ├── phenotype/phenotype_formatting.py
│       └── ...
```

---

## Task 1 — MWE-Level Dry-Run Tests (Start Here)

Validate that the route3 pipeline passes all 21 dry-run tests before making structural changes.

**Current state**: The 21 scripts in `dryrun/mwe/NN_<step>.sh` reference `pipeline/` notebooks. They must be updated to call `route3/notebooks/` instead, then run and compared against expected output.

### Steps

**1. Update all MWE scripts** — for each `dryrun/mwe/NN_<step>.sh`, replace the `PIPE` variable:

```bash
# Old:
PIPE=/home/user/xqtl-protocol/pipeline
# New:
PIPE=/home/user/xqtl-protocol/renovated_code/snakemake/route3/notebooks
```

The notebook filename is in the script's header comment, e.g.:
`# MWE: RNA_calling.ipynb :: fastqc` → file is `$PIPE/RNA_calling.ipynb`

**2. Run all scripts and compare against expected output**:

```bash
cd /home/user/xqtl-protocol
for sh in renovated_code/snakemake/dryrun/mwe/*.sh; do
    num=$(basename "$sh" | cut -c1-2)
    step=$(basename "$sh" .sh | cut -c4-)
    bash "$sh" > /tmp/mwe_${num}.txt 2>&1
    expected="renovated_code/snakemake/dryrun/output/${num}_${step}.txt"
    if diff -q "$expected" /tmp/mwe_${num}.txt > /dev/null 2>&1; then
        echo "PASS ${num}_${step}"
    else
        echo "FAIL ${num}_${step}"
        diff "$expected" /tmp/mwe_${num}.txt
    fi
done
```

**3. Fix any failures.** Common causes:

- Route3 notebook `[global]` block missing `renovated_code_dir` parameter — fix in `transform_notebooks.py`, then regenerate:
  ```bash
  python3 renovated_code/snakemake/route3/transform_notebooks.py
  ```
- Script runner set to `bash` when it should be `Rscript` or `python3` — fix the `make_bash_block()` call for that step in `transform_notebooks.py`
- SoS step parameters in the route3 notebook not matching what the MWE script passes — align parameters in the transform function

**Acceptance criteria**: all 21 scripts exit 0 with empty diffs (minor whitespace differences are acceptable).

---

## Task 2 — Strict 1-1 Mini-Protocol ↔ Snakemake Mapping + MemVerge Parallelization

### What is a mini-protocol?

A mini-protocol is a top-level Jupyter notebook (e.g., `code/molecular_phenotypes/bulk_expression.ipynb`) with numbered steps (i, ii, iii…). Each step has a companion code cell running one `sos run` command for one SoS workflow step. These are the **ground truth** for what the pipeline should do.

Example — `bulk_expression.ipynb` maps to `01_molecular_phenotypes.smk`:

| Step | Mini-protocol header | SoS step name | Expected Snakemake rule |
|------|---------------------|---------------|------------------------|
| i | FastQC quality summary | `fastqc` | `fastqc` |
| ii | Cut adaptor (Optional) | `fastp_trim_adaptor` | `fastp_trim_adaptor` |
| iii | STAR alignment + Picard QC | `STAR_align` | `star_align` |
| iv | Gene-level RNA expression via rnaseqc | `rnaseqc_call` | `rnaseqc_call` |
| v | Transcript-level expression via RSEM | `rsem_call` | `rsem_call` |
| vi | Multi-sample RNA-seq QC | `qc` | `bulk_expression_qc` |
| vii | Read count normalization | `normalize` | `bulk_expression_normalization` |

**The 1-1 requirement operates at two levels**:

1. **File level**: each `.smk` module file corresponds to exactly one mini-protocol notebook
2. **Rule level**: each numbered step in the mini-protocol corresponds to exactly one Snakemake rule — no batching, no skipping

### Current file-level mapping to audit and enforce

| SMK module | Mini-protocol notebook |
|---|---|
| `01_molecular_phenotypes.smk` | `code/molecular_phenotypes/bulk_expression.ipynb` |
| `02_genotype_preprocessing.smk` | `code/data_preprocessing/genotype_preprocessing.ipynb` |
| `03_sample_qc_pca.smk` | `code/data_preprocessing/genotype_preprocessing.ipynb` (PCA section) |
| `04_phenotype_covariate_prep.smk` | `code/data_preprocessing/phenotype_preprocessing.ipynb` + `covariate_preprocessing.ipynb` |
| `05_association_testing.smk` | `code/association_scan/qtl_association_testing.ipynb` |
| `06_univariate_finemapping.smk` | `code/mnm_analysis/mnm_miniprotocol.ipynb` |

Note: `03` and `04` currently span multiple mini-protocols. Audit whether these should be split into separate `.smk` files to achieve clean 1-1 at the file level.

### What is MemVerge-style parallelization?

Currently, one Snakemake rule submits one `sos run` call that **internally loops** over units (regions, chromosomes, phenotype groups). This hides parallelism from Snakemake — if one region fails, the entire job must restart from the beginning.

In the target design, **Snakemake owns the parallelism**: each atomic unit becomes its own wildcard-driven rule instance, submitted as a separate cluster job. Individual failures are surgically restartable.

Real example from production usage — what the target looks like:

```bash
# Each TADB region = its own independently restartable cluster job
sos run pipeline/mnm_regression.ipynb fsusie --region_name TADB_1003 --name Astro ...
sos run pipeline/mnm_regression.ipynb fsusie --region_name TADB_1008 --name Astro ...
sos run pipeline/mnm_regression.ipynb fsusie --region_name TADB_1032 --name Astro ...
```

### Implementation pattern (checkpoint + per-unit rule)

```python
# Step 1: Checkpoint reads the region list and creates one flag file per region.
# Snakemake re-evaluates the DAG after this runs to discover all regions.
checkpoint enumerate_regions:
    input:
        region_list = "{cwd}/finemapping/{theme}/region_list.txt"
    output:
        directory("{cwd}/finemapping/{theme}/regions/")
    shell:
        "mkdir -p {output}; "
        "awk 'NR>1{{print $1}}' {input} | xargs -I{{}} touch {output}/{{}}.flag"

# Step 2: Aggregator function — collects all per-region done files.
# Called in the 'all' rule or a downstream rule's input.
def all_region_outputs(wildcards):
    checkpoints.enumerate_regions.get(**wildcards)
    region_dir = Path(f"{wildcards.cwd}/finemapping/{wildcards.theme}/regions/")
    return expand(
        "{cwd}/finemapping/{theme}/susie_twas/{region}.done",
        region=[f.stem for f in region_dir.glob("*.flag")],
        **wildcards,
    )

# Step 3: Per-unit rule — one cluster job submitted per region.
rule susie_twas_region:
    input:
        flag   = "{cwd}/finemapping/{theme}/regions/{region}.flag",
        geno   = "{cwd}/data_preprocessing/genotype/xqtl_protocol_data.plink_qc.genotype_by_chrom_files.txt",
        pheno  = "{cwd}/data_preprocessing/{theme}/phenotype_data/{theme}.phenotype_by_chrom_files.txt",
        hidden_factors = lambda wc: get_hidden_factors(wc),
    output:
        done = "{cwd}/finemapping/{theme}/susie_twas/{region}.done"
    params:
        notebooks_dir = ROUTE3_NOTEBOOKS,
        renovated_dir = RENOVATED_CODE,
        container     = config["containers"]["susie"],
        outdir        = "{cwd}/finemapping/{theme}/susie_twas",
        # ... other params from config
    threads: config["resources"]["finemapping"]["threads"]
    resources:
        mem_mb  = config["resources"]["finemapping"]["mem_mb"],
        runtime = config["resources"]["finemapping"]["runtime"],
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.region}
        sos run {params.notebooks_dir}/mnm_regression.ipynb susie_twas \
            --region_name {wildcards.region} \
            --name {wildcards.theme} \
            --genoFile {input.geno} \
            --phenoFile {input.pheno} \
            --covFile {input.hidden_factors} \
            --cwd {params.outdir}/{wildcards.region} \
            --container {params.container} \
            --renovated-code-dir {params.renovated_dir} \
            --numThreads {threads}
        touch {output.done}
        """
```

### Unit of parallelism per step

Read each mini-protocol MWE command to confirm the natural unit. Current best mapping:

| Step category | Unit | Wildcard |
|---|---|---|
| Finemapping (`susie_twas`, `fsusie`) | TADB region ID | `{region}` |
| QTL association (`tensorqtl_cis`) | Chromosome | `{chrom}` |
| Phenotype splitting (`phenotype_by_chrom`) | Chromosome | `{chrom}` |
| Per-tissue steps (QC, normalization, hidden factors) | Tissue/theme | `{theme}` (already wildcarded) |

### Audit steps

**1. Find rules that still use SoS-internal looping:**

```bash
grep -n "region_list\|pheno_list\|chrom_list\|phenoFile.*list" \
    renovated_code/snakemake/rules/sos/*.smk
```

**2. Read each mini-protocol in `code/` and produce the full audit table:**

| Mini-protocol | Step | SoS step | Snakemake rule | Parallelism unit | Status |
|---|---|---|---|---|---|
| `bulk_expression.ipynb` | i | `fastqc` | `fastqc` | `{theme}` | ✅ |
| `bulk_expression.ipynb` | iii | `STAR_align` | `star_align` | `{theme}` | ✅ or ⚠️ |
| `mnm_miniprotocol.ipynb` | fsusie | `susie_twas` | `susie_twas` | `{region}` | ⚠️ SoS-internal |
| ... | | | | | |

Status values: `✅` 1-1 + Snakemake-parallel | `⚠️` 1-1 but SoS-internal | `❌` missing rule

**3. Refactor all ⚠️ rules** using the checkpoint pattern above.

**4. Validate:** `snakemake -n` must report N jobs = N regions/chroms for multi-unit steps, not 1.

---

## Task 3 — Extract All Inline Blocks to Modular Scripts

### Background

Route3 notebooks are generated by `transform_notebooks.py`. It replaces inline `R:` and `python:` task blocks with calls to external scripts in `renovated_code/`. Some steps are not yet extracted and still contain inline code in the route3 notebooks.

### Find remaining inline blocks

Run this script from the repo root:

```bash
python3 - << 'EOF'
import json, glob, re
from pathlib import Path

NOTEBOOKS = sorted(glob.glob(
    "renovated_code/snakemake/route3/notebooks/*.ipynb"
))

# Collect steps currently called by active Snakemake rules
ACTIVE_STEPS = set()
for smk in glob.glob("renovated_code/snakemake/rules/sos/*.smk"):
    with open(smk) as f:
        content = f.read()
    for m in re.finditer(r"sos run .*?\.ipynb\s+(\S+)", content):
        ACTIVE_STEPS.add(m.group(1))

print(f"Active Snakemake steps: {sorted(ACTIVE_STEPS)}\n")

for nb_path in NOTEBOOKS:
    nb_name = Path(nb_path).name
    with open(nb_path) as f:
        nb = json.load(f)
    current_step = None
    for cell in nb["cells"]:
        src = "".join(cell.get("source", []))
        if cell["cell_type"] == "code":
            step_m = re.match(r"^\[(\S+)\]", src)
            if step_m:
                current_step = step_m.group(1)
            if re.search(r"^(R|python|bash):\s*$", src, re.MULTILINE):
                priority = "HIGH" if current_step in ACTIVE_STEPS else "low"
                first_line = next(
                    (ln for ln in src.splitlines() if re.match(r"(R|python|bash):", ln)),
                    ""
                )
                print(f"[{priority}]  {nb_name}  step={current_step}  block_type={first_line!r}")
EOF
```

Focus on **HIGH** priority blocks first — those are called by active Snakemake rules and block production use.

### Extraction pattern for task blocks

For each inline block in step `[step_name]`:

1. **Find or create the target script** in `renovated_code/<category>/<tool>.{R,py}`:
   - If the script exists, add a new `--step step_name` dispatch branch to its argument parser
   - If not, create it with `argparse` (Python) or `optparse` (R) boilerplate

2. **Add or update the transform in `transform_notebooks.py`**:

   ```python
   # For an R step:
   make_bash_block(
       step_name="step_name",
       runner="rscript",
       script_rel_path="category/tool.R",
       params=[
           "--param1 {param1}",
           "--param2 {param2}",
           "--numThreads {numThreads}",
       ],
   )

   # For a Python step:
   make_bash_block(
       step_name="step_name",
       runner="python3",
       script_rel_path="category/tool.py",
       params=[
           "--param1 {param1}",
           "--numThreads {numThreads}",
       ],
   )
   ```

   For R and Python runners, SoS handles containers via `container = container` in the task block header — do **not** add a `--container` flag.

3. **Regenerate route3 notebooks:**
   ```bash
   python3 renovated_code/snakemake/route3/transform_notebooks.py
   ```

4. **Re-run the relevant MWE test** to confirm the step still passes:
   ```bash
   bash renovated_code/snakemake/dryrun/mwe/NN_<step>.sh
   ```

### Global session extraction (complex notebooks only)

For notebooks whose `[global]` section contains **non-trivial logic** — file discovery loops, conditional branches, dynamic list construction — extract it to a resolver script. Do **not** apply this to notebooks with simple global sections (just imports and parameter defaults).

1. Extract global logic to `renovated_code/<category>/<notebook>_params.py`
2. The script outputs a `params.json`:
   ```json
   {"phenoFile": ["/path/chr1.bed", "/path/chr2.bed"], "covFile": [...]}
   ```
3. The Snakemake rule reads the JSON:
   ```python
   rule my_step:
       input:
           params_json = "{cwd}/params/{theme}_params.json",
           ...
       params:
           p = lambda wc, input: json.load(open(input.params_json))
       shell:
           "sos run {params.notebooks_dir}/notebook.ipynb step "
           "--phenoFile {params.p[phenoFile]} ..."
   ```

### Simple bash wrappers — no changes needed

Scripts like `RNA_calling.sh`, `VCF_QC.sh`, `genotype_formatting.sh` wrap CLI tools (STAR, BCFtools, PLINK) directly. Leave them unchanged.

### Goal

Every extracted script must be independently runnable:

```bash
Rscript renovated_code/category/tool.R --help
python3 renovated_code/category/tool.py --help
```

This enables `pytest` / `testthat` unit tests in the future without requiring the full pipeline environment, containers, or real data.

---

## Validation Commands

These are the raw commands to use for linting and checking — no environment-specific slash commands needed:

```bash
# Lint the Snakemake pipeline (catches syntax/rule errors)
snakemake --snakefile renovated_code/snakemake/Snakefile \
          --configfile renovated_code/snakemake/config.yaml \
          --lint

# Dry-run the route3 pipeline to check the DAG
snakemake --snakefile renovated_code/snakemake/route3/Snakefile \
          --configfile renovated_code/snakemake/config.yaml \
          -n --quiet

# Find stale pipeline/ references in rule files
grep -rn "pipeline/" renovated_code/snakemake/rules/ --include="*.smk" | grep -v "^\s*#"

# Scan for hardcoded absolute paths outside of config expansions
grep -rn "^[^#]*/home/\|^[^#]*/tmp/\|^[^#]*/opt/" \
    renovated_code/snakemake/rules/ --include="*.smk"

# Regenerate all route3 notebooks
python3 renovated_code/snakemake/route3/transform_notebooks.py
```

---

## Commit Conventions

Format: `<type>(<scope>): <summary>`

| Scope | When to use |
|-------|-------------|
| `route3` | `transform_notebooks.py` or `route3/notebooks/` |
| `smk` | Any `.smk` or `Snakefile` change |
| `scripts` | `renovated_code/` R/Python/bash scripts |
| `dryrun` | `dryrun/mwe/` or `dryrun/output/` |

Always push to the current feature branch. Never push to `main` or `master`.
