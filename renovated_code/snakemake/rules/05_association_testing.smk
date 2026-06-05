# ============================================================
# Rule Module 05: QTL Association Testing (TensorQTL)
# ============================================================
# Covers: CIS-QTL association testing with TensorQTL
#
# Mirrors: eQTL_analysis_commands.ipynb stage: TensorQTL
#
# SoS notebooks called:
#   - pipeline/TensorQTL.ipynb (cis)
#
# Note: get_hidden_factors() and get_phenotype_bed() are defined in the main
# Snakefile and are available to all included rule modules.
#
# TensorQTL `cis` workflow runs two internal sub-steps automatically:
#   cis_1: per-chromosome nominal scan + permutation → parquet + TSV files
#   cis_2: aggregate per-chromosome results, compute q-values → significance TSV
# There is NO separate cis_postprocess step to call.
# ============================================================

# ------------------------------------
# Step 5.1 — TensorQTL CIS association scan
# ------------------------------------
# Inputs:
#   --genotype-file: path to genotype_by_chrom_files.txt (list of per-chrom plinks)
#   --phenotype-file: path to {theme}.phenotype_by_chrom_files.txt (list of per-chrom BEDs)
#   --covariate-file: hidden factor file (Marchenko PCA or PEER)
#
# Output: uses a sentinel done file because TensorQTL produces many files
# whose names depend on the input file basenames. The done file is written
# only after cis_2 completes (both nominal + significance aggregation).
rule tensorqtl_cis:
    """Run TensorQTL cis-QTL nominal + permutation scan across all chromosomes."""
    input:
        geno_list      = lambda wc: get_genotype_chrom_list(),
        phenotype_bed  = lambda wc: get_phenotype_bed(wc),
        pheno_list     = lambda wc: get_phenotype_chrom_list(wc.theme),
        hidden_factors = lambda wc: get_hidden_factors(wc),
    output:
        done = "{cwd}/association_scan/{theme}/TensorQTL/.done_cis",
    params:
        tensorqtl_notebook = get_pipeline_notebook_path("TensorQTL.ipynb"),
        compat_python  = COMPAT_PYTHON,
        outdir       = "{cwd}/association_scan/{theme}/TensorQTL",
        normalized_geno_manifest = lambda wc: (
            f"{config['cwd']}/association_scan/{wc.theme}/TensorQTL/"
            f"{PLINK_QC_BASENAME}.genotype_by_chrom_files.normalized.txt"
        ),
        normalized_pheno_manifest = lambda wc: (
            f"{config['cwd']}/association_scan/{wc.theme}/TensorQTL/"
            f"{get_phenotype_base(wc.theme)}.phenotype_by_chrom_files.normalized.txt"
        ),
        cis_window   = config["association"]["cis_window"],
        mac          = config["association"]["mac_threshold"],
        maf          = config["association"]["maf_threshold"],
        dry_run      = DRY_RUN_SOS,
    threads: config["resources"]["tensorqtl"]["threads"]
    resources:
        mem_mb   = config["resources"]["tensorqtl"]["mem_mb"],
        runtime = config["resources"]["tensorqtl"]["runtime"],
    shell:
        """
        mkdir -p {params.outdir}
        awk 'BEGIN{{FS=OFS="\\t"}} NR==1 {{print "#id","#path"; next}} {{$1=$1; sub(/^chr/, "", $1); print $1, $2}}' {input.geno_list} > {params.normalized_geno_manifest}
        awk 'BEGIN{{FS=OFS="\\t"}} NR==1 {{print "#id","#path"; next}} {{$1=$1; sub(/^chr/, "", $1); print $1, $2}}' {input.pheno_list} > {params.normalized_pheno_manifest}
        geno_rows="$(awk 'END {{print NR-1}}' {params.normalized_geno_manifest})"
        pheno_rows="$(awk 'END {{print NR-1}}' {params.normalized_pheno_manifest})"
        chrom_args="$(awk 'NR>1 {{print $1}}' {params.normalized_geno_manifest} | paste -sd ' ' -)"
        if [ -z "$chrom_args" ]; then
            echo "ERROR: normalized genotype manifest is empty" >&2
            exit 1
        fi
        if [ "$geno_rows" -eq 1 ] && [ "$pheno_rows" -eq 1 ]; then
            geno_bed="$(awk 'NR==2 {{print $2}}' {params.normalized_geno_manifest})"
            XQTL_PATCH_TENSORQTL_SORT=1 PYTHONPATH="{params.compat_python}${{PYTHONPATH:+:$PYTHONPATH}}" \
            sos run {params.tensorqtl_notebook} cis {params.dry_run} \
                --cwd {params.outdir} \
                --genotype-file "$geno_bed" \
                --phenotype-file {input.phenotype_bed} \
                --chromosome $chrom_args \
                --covariate-file {input.hidden_factors} \
                --window {params.cis_window} \
                --MAC {params.mac} \
                --maf-threshold {params.maf} \
                --numThreads {threads}
        else
            XQTL_PATCH_TENSORQTL_SORT=1 PYTHONPATH="{params.compat_python}${{PYTHONPATH:+:$PYTHONPATH}}" \
            sos run {params.tensorqtl_notebook} cis {params.dry_run} \
                --cwd {params.outdir} \
                --genotype-file {params.normalized_geno_manifest} \
                --phenotype-file {params.normalized_pheno_manifest} \
                --chromosome $chrom_args \
                --covariate-file {input.hidden_factors} \
                --window {params.cis_window} \
                --MAC {params.mac} \
                --maf-threshold {params.maf} \
                --numThreads {threads}
        fi
        touch {output.done}
        """
