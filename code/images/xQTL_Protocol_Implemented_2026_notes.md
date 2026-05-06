# Updated FunGen-xQTL workflow plot notes

## Relationship to the existing Lucidchart figure

`xQTL_Protocol_Implemented_2026.svg` is now a wholesale visual extension of the existing protocol figure. The generator embeds the current Lucidchart-derived PNG, `xQTL_Protocol_Feb_2024.png`, directly into the SVG as a base64 image layer and appends a right-side panel for components that should be added or made more explicit. This makes the SVG self-contained, so opening it outside the repository still shows the original figure plus the added panel.

The Lucidchart JSON export, `Protocol _2024_Mar.json`, is still used as a source audit. The script checks for expected labels from the original schematic before writing the extended SVG, so the generated asset remains tied to the existing Lucidchart vocabulary.

## What was added to the existing plot

The appended panel adds repository components that are missing from, or not explicit enough in, the existing visual schematic:

- SuSiE-enloc enrichment for xQTL-GWAS enrichment priors before pairwise colocalization.
- cTWAS workflow after TWAS and variant selection for candidate genes.
- ColocBoost as an additional GWAS/xQTL colocalization method.
- EOO annotation enrichment with chromosome jackknife validation.
- Pathway / GSEA enrichment for interpreting prioritized genes.
- GREGOR and sLDSC enrichment workflows.
- xQTL modifier score / EMS training and prediction.

## Why this version is better for the requested purpose

This version keeps the existing Lucidchart figure visually intact rather than replacing it with a simplified redesign. It is better for the requested purpose because reviewers can still recognize the original protocol schematic while seeing a clearly marked add-on panel with the missing implemented repository components.

## Limitations

Because the Lucidchart JSON export in this repository does not include shape coordinates or style geometry, the generator cannot reconstruct the original vector drawing from JSON alone. It therefore preserves the existing PNG as the base figure and uses the JSON only to audit expected source labels.
