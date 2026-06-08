# xQTL Modular SoS Handoff

## Human Note

This branch replaces the old inline implementation notebooks with modular SoS notebooks in the canonical `pipeline/` layer. The separate modular notebook staging directories were removed, so reviewers should see one notebook source of truth rather than duplicated notebook copies.

## Current Layout

- `pipeline/*.ipynb` are symlinks into the corresponding implementation notebooks under `code/`.
- The promoted notebooks call modular scripts under `renovated_code/script/`.
- `renovated_code/snakemake/modular_sos/Snakefile` runs `sos run pipeline/<notebook>.ipynb <step>`.
- `renovated_code/snakemake/Snakefile` also defaults modular SoS notebook redirection to `pipeline/`.
- Removed duplicate notebook staging layers:
  - `renovated_code/snakemake/modular_sos/notebooks/`
  - `renovated_code/notebook/modular_sos/`
- Removed the generated-notebook transformer. Notebook changes now live directly in the canonical notebook layer.

## Runtime Bundle

The minimal runtime bundle for MWE and MESA remains:

- `renovated_code/snakemake/modular_sos/Snakefile`
- `renovated_code/snakemake/modular_sos/rules/00_*.smk` through `06_*.smk`
- canonical notebooks reached via `pipeline/`
- modular scripts under `renovated_code/script/`
- tests under `renovated_code/snakemake/modular_sos/tests/`

## Validation Commands

Static checks used for this cleanup:

```bash
jq empty pipeline/GWAS_QC.ipynb pipeline/PCA.ipynb pipeline/RNA_calling.ipynb \
  pipeline/TensorQTL.ipynb pipeline/VCF_QC.ipynb pipeline/bulk_expression_QC.ipynb \
  pipeline/bulk_expression_normalization.ipynb pipeline/covariate_formatting.ipynb \
  pipeline/covariate_hidden_factor.ipynb pipeline/genotype_formatting.ipynb \
  pipeline/mnm_regression.ipynb pipeline/phenotype_formatting.ipynb \
  pipeline/phenotype_imputation.ipynb pipeline/rss_analysis.ipynb

git diff --check
```

Runtime checks to rerun before merge:

```bash
renovated_code/snakemake/modular_sos/tests/run_mwe_xqtl_core.sh \
  --mwe-data ../mwe_data \
  --run-tag modular_sos_mwe_core \
  --cores 1

renovated_code/snakemake/modular_sos/tests/run_nontrivial_tensorqtl_susie.sh \
  --run-tag modular_sos_nontrivial_compare \
  --num-threads 4
```

For MESA or other downstream data:

```bash
snakemake \
  --snakefile renovated_code/snakemake/modular_sos/Snakefile \
  --configfile path/to/modular_sos.config.yaml \
  --cores 1 \
  xqtl_core
```

Use target `all` only when fine-mapping plot generation is required.

## Fidelity

Fidelity verdict: `faithful`. The Snakemake orchestration remains notebook-first, the notebook stage graph is preserved, and implementation-heavy task bodies now delegate to modular scripts instead of maintaining duplicated inline code.
