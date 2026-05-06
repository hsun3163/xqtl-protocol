# Updated FunGen-xQTL workflow plot notes

## Relationship to the existing Lucidchart figure

The generated `xQTL_Protocol_Implemented_2026.svg` is not intended to replace the detailed Lucidchart schematic exported as `Protocol _2024_Mar.json`. The Lucidchart export is the richer source diagram: it contains multiple pages and dozens of labels covering fine-grained protocol concepts, including data sources, omics contexts, QC, association testing, fine-mapping, integration, prioritization, and Alzheimer’s disease follow-up.

The generated SVG is a condensed implementation view. It groups related Lucidchart nodes into fewer boxes so readers can quickly see which broad workflow families are represented in the current repository and how downstream integrative analyses follow association testing and fine-mapping.

## What changed relative to the existing plot

- The detailed Lucidchart boxes for existing GWAS/xQTL summary statistics, reference genotype data, LD reference panels, and context-specific xQTL assays are grouped into an `Inputs & reference data` panel.
- Genotype QC, PCA, kinship, allele/build harmonization, phenotype QC, imputation, and covariate inference are grouped into a `QC, harmonization & covariates` panel.
- xQTL association, quantile QTL, univariate fine-mapping/TWAS weights, multivariate or multi-gene fine-mapping, and omics prediction databases are grouped into a `QTL association & fine-mapping` panel.
- The downstream integrative-analysis portion is made more explicit with a `GWAS / trait integration` panel covering SuSiE-enloc enrichment, pairwise colocalization, TWAS to cTWAS to MR, and ColocBoost / coloc post-processing.
- Validation and prioritization are separated into a bottom panel covering annotation/pathway/heredity enrichment, EMS, prioritization outputs, and Alzheimer’s disease follow-up hypotheses.

## Why this can be useful

The existing Lucidchart figure is better for a full conceptual overview. The generated SVG is better only for a narrower purpose: it is smaller, reproducible from source control, easier to edit programmatically, and emphasizes implemented downstream analysis families after QTL association and fine-mapping.

## Limitations

This generated SVG deliberately omits many fine-grained Lucidchart labels and citations to keep the figure readable. It should therefore be treated as a documentation companion or implementation overlay, not as the authoritative protocol schematic.
