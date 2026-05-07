# Updated FunGen-xQTL workflow plot notes

## Relationship to the existing Lucidchart figure

`xQTL_Protocol_Implemented_2026.tex` is a TikZ implementation of the updated workflow figure. It preserves the existing Lucidchart-derived protocol PNG, `xQTL_Protocol_Feb_2024.png`, as the base layer and draws the missing or newly explicit repository components in TikZ as a right-side add-on panel.

The Lucidchart JSON export, `Protocol _2024_Mar.json`, remains the source vocabulary for the original schematic. The TikZ file does not try to redraw the original schematic from JSON because this JSON export does not include usable shape coordinates or style geometry; instead, it keeps the original visual intact and implements the requested additions as native TikZ nodes and arrows.

## What was added to the existing plot

The TikZ panel adds repository components that are missing from, or not explicit enough in, the existing visual schematic:

- SuSiE-enloc enrichment for xQTL-GWAS enrichment priors before pairwise colocalization.
- cTWAS workflow after TWAS and variant selection for candidate genes.
- ColocBoost as an additional GWAS/xQTL colocalization method.
- EOO annotation enrichment with chromosome jackknife validation.
- Pathway / GSEA enrichment for interpreting prioritized genes.
- GREGOR and sLDSC enrichment workflows.
- xQTL modifier score / EMS training and prediction.

## How to compile

From the repository root, run:

```bash
pdflatex code/images/xQTL_Protocol_Implemented_2026.tex
```

Alternatively, compile the file from `code/images/`; the `\graphicspath` directive supports both locations.
