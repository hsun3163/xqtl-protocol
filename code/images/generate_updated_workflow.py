"""Generate a condensed implementation-view workflow SVG.

This figure is a companion to the detailed Lucidchart protocol schematic in
``Protocol _2024_Mar.json``. It intentionally groups several Lucidchart nodes
into broader repository-implementation categories so the downstream workflow is
easier to read in documentation and presentations.
"""

from pathlib import Path
import html
import json
import re

OUT = Path(__file__).with_name('xQTL_Protocol_Implemented_2026.svg')
LUCIDCHART_SOURCE = Path(__file__).with_name('Protocol _2024_Mar.json')
W, H = 1640, 1080

EXPECTED_LUCIDCHART_LABELS = [
    'Existing GWAS & xQTL Summary Statistics',
    'Context Specfic xQTL',
    'Genotype Quality Control',
    'Omics-data Quantification',
    'Genome-wide xQTL Association',
    'Genome-wide xQTL Fine-mapping',
    'Colocalization',
    'TWAS',
    'Causal TWAS',
    'Mendelian Randomization',
    'Multi-Context xQTL',
    'Variant-level Prioritization',
    "Alzheimer's Disease",
]

COLORS = {
    'bg': '#fbfbf7', 'panel': '#ffffff', 'stroke': '#344054', 'muted': '#667085',
    'input': '#dbeafe', 'process': '#dcfce7', 'integrate': '#fee2e2',
    'validate': '#fef3c7', 'output': '#e0e7ff', 'line': '#344054', 'dash': '#94a3b8'
}


def lucidchart_labels():
    if not LUCIDCHART_SOURCE.exists():
        raise FileNotFoundError(f'Missing Lucidchart source: {LUCIDCHART_SOURCE}')
    data = json.loads(LUCIDCHART_SOURCE.read_text())
    labels = set()
    for page in data.get('pages', []):
        for shape in page.get('items', {}).get('shapes', []):
            text = ' '.join(area.get('text', '') for area in shape.get('textAreas', []))
            text = re.sub(r'\s+', ' ', text).strip()
            if text:
                labels.add(text)
    return labels


def assert_lucidchart_source_coverage():
    labels = lucidchart_labels()
    missing = [label for label in EXPECTED_LUCIDCHART_LABELS if label not in labels]
    if missing:
        raise ValueError('Expected Lucidchart labels are missing: ' + ', '.join(missing))
    print(f'Audited {len(labels)} Lucidchart labels from {LUCIDCHART_SOURCE.name}')


def wrap(text, max_chars):
    words = text.split()
    lines, cur = [], ''
    for word in words:
        if len(cur) + len(word) + (1 if cur else 0) <= max_chars:
            cur = f"{cur} {word}".strip()
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


def text_block(x, y, lines, size=18, weight='600', fill='#111827', anchor='middle', line_height=1.2):
    out = []
    for idx, line in enumerate(lines):
        dy = idx * size * line_height
        out.append(f'<text x="{x}" y="{y + dy:.1f}" text-anchor="{anchor}" font-size="{size}" font-weight="{weight}" fill="{fill}">{html.escape(line)}</text>')
    return '\n'.join(out)


def box(x, y, w, h, title, subtitle='', fill='process', stroke=None, radius=14, title_size=19, sub_size=15, max_chars=24):
    stroke = stroke or COLORS['stroke']
    title_lines = wrap(title, max_chars)
    sub_lines = wrap(subtitle, max_chars + 8) if subtitle else []
    line_count = len(title_lines) + len(sub_lines)
    start_y = y + h / 2 - ((line_count - 1) * title_size * 1.15) / 2 + title_size / 3
    parts = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{COLORS[fill]}" stroke="{stroke}" stroke-width="2"/>']
    parts.append(text_block(x + w / 2, start_y, title_lines, size=title_size, weight='700'))
    if sub_lines:
        parts.append(text_block(x + w / 2, start_y + len(title_lines) * title_size * 1.18 + 6, sub_lines, size=sub_size, weight='500', fill=COLORS['muted']))
    return '\n'.join(parts)


def arrow(x1, y1, x2, y2, dash=False):
    style = 'stroke-dasharray="8 8"' if dash else ''
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{COLORS["line"]}" stroke-width="3" marker-end="url(#arrow)" {style}/>'


def panel(x, y, w, h, title):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="{COLORS["panel"]}" stroke="{COLORS["dash"]}" stroke-width="2.5" stroke-dasharray="14 10"/>\n' + text_block(x + w / 2, y + 34, [title], size=22, weight='800', fill='#1f2937')


svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
    '<defs><marker id="arrow" markerWidth="12" markerHeight="12" refX="10" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,6 L10,3 z" fill="#344054"/></marker></defs>',
    f'<rect width="100%" height="100%" fill="{COLORS["bg"]}"/>',
]
svg.append(text_block(W / 2, 48, ['FunGen-xQTL Protocol Workflow: implemented 2026 repository view'], size=30, weight='900'))
svg.append(text_block(W / 2, 80, ['Condensed implementation view audited against Protocol _2024_Mar.json; not a replacement for the detailed schematic.'], size=16, weight='500', fill=COLORS['muted']))

svg.append(panel(50, 115, 340, 425, 'Inputs & reference data'))
svg.append(panel(430, 115, 340, 425, 'QC, harmonization & covariates'))
svg.append(panel(810, 115, 360, 425, 'QTL association & fine-mapping'))
svg.append(panel(1210, 115, 380, 425, 'GWAS / trait integration'))
svg.append(panel(430, 600, 1160, 360, 'Validation, prioritization & scores'))

svg.append(box(90, 165, 260, 72, 'GWAS meta-analysis summary statistics', fill='input'))
svg.append(box(90, 258, 260, 72, 'xQTL summary statistics', fill='input'))
svg.append(box(90, 351, 260, 72, 'Reference genotype and LD panels', fill='input'))
svg.append(box(90, 444, 260, 72, 'Context-specific omics', 'expression, splicing, methylation, APA, chromatin, protein / metabolites', fill='input', max_chars=26))

svg.append(box(470, 158, 260, 70, 'Genotype QC, PCA and kinship', fill='process'))
svg.append(box(470, 252, 260, 80, 'Allele/build harmonization and Z-score correction', fill='process'))
svg.append(box(470, 356, 260, 80, 'Phenotype QC, normalization and imputation', fill='process'))
svg.append(box(470, 460, 260, 58, 'Known + hidden covariates', fill='process'))

svg.append(box(850, 150, 280, 72, 'xQTL association testing', 'TensorQTL; quantile QTL', fill='validate'))
svg.append(box(850, 246, 280, 82, 'Univariate fine-mapping and TWAS weights', 'SuSiE, RSS, fSuSiE', fill='validate'))
svg.append(box(850, 352, 280, 82, 'Multivariate / multi-gene fine-mapping', 'mvSuSiE, mr.mash, TAD windows', fill='validate'))
svg.append(box(850, 458, 280, 58, 'Omics prediction databases', fill='output'))

svg.append(box(1250, 150, 300, 72, 'SuSiE-enloc enrichment', 'global xQTL-GWAS priors', fill='integrate'))
svg.append(box(1250, 246, 300, 72, 'Pairwise colocalization', 'shared causal variants', fill='integrate'))
svg.append(box(1250, 342, 300, 82, 'TWAS → cTWAS → MR', 'gene-trait association and causal follow-up', fill='integrate'))
svg.append(box(1250, 458, 300, 58, 'ColocBoost and coloc post-processing', fill='integrate'))

svg.append(box(475, 665, 235, 72, 'EOO annotation enrichment', 'LOCO jackknife', fill='process'))
svg.append(box(745, 665, 235, 72, 'Pathway / gene set enrichment', 'GSEA multi-database', fill='process'))
svg.append(box(1015, 665, 235, 72, 'GREGOR and sLDSC', 'variant and heritability enrichment', fill='process'))
svg.append(box(1285, 665, 235, 72, 'xQTL modifier score', 'EMS training and prediction', fill='process'))
svg.append(box(600, 810, 260, 76, 'Variant, locus, gene and region prioritization', fill='output'))
svg.append(box(1015, 810, 260, 76, "Alzheimer's disease biology and follow-up hypotheses", fill='output'))

for y in [201, 294, 390, 480]:
    svg.append(arrow(350, y, 470, y))
for y1, y2 in [(193, 186), (292, 286), (396, 393), (489, 487)]:
    svg.append(arrow(730, y1, 850, y2))
for y1, y2 in [(186, 186), (287, 282), (392, 383), (487, 487)]:
    svg.append(arrow(1130, y1, 1250, y2))
for y1, y2 in [(222, 246), (328, 352), (434, 458)]:
    svg.append(arrow(980, y1, 980, y2))
for y1, y2 in [(222, 246), (318, 342), (424, 458)]:
    svg.append(arrow(1400, y1, 1400, y2))
for x in [1325, 1400, 1475, 980, 900, 1060]:
    svg.append(arrow(x, 540, x, 640, dash=True))
svg.append(arrow(592, 737, 706, 810))
svg.append(arrow(862, 737, 770, 810))
svg.append(arrow(1132, 737, 1145, 810))
svg.append(arrow(1402, 737, 1145, 810))
svg.append(arrow(860, 848, 1015, 848))

legend_y = 1000
items = [
    ('input', 'Inputs'),
    ('process', 'Implemented processing / validation'),
    ('validate', 'Association and fine-mapping'),
    ('integrate', 'Implemented GWAS integration'),
    ('output', 'Outputs'),
]
x = 80
for key, label in items:
    svg.append(f'<rect x="{x}" y="{legend_y - 18}" width="28" height="18" rx="4" fill="{COLORS[key]}" stroke="{COLORS["stroke"]}"/>')
    svg.append(f'<text x="{x + 38}" y="{legend_y - 4}" font-size="15" font-weight="600" fill="#374151">{html.escape(label)}</text>')
    x += 300 if key != 'input' else 180

svg.append(text_block(W / 2, 1045, [f'Source audit: {LUCIDCHART_SOURCE.name}; grouped for readability and implementation status.'], size=14, weight='500', fill=COLORS['muted']))
svg.append('</svg>')
assert_lucidchart_source_coverage()
OUT.write_text('\n'.join(svg))
print(OUT)
