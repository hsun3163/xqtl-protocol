"""Generate a Lucidchart-based workflow SVG with repository additions.

The generated figure keeps the existing protocol figure intact by embedding the
current Lucidchart-derived PNG (`xQTL_Protocol_Feb_2024.png`) as the base layer.
It then adds a right-side panel for repository components that are not explicit
in the existing Lucidchart export (`Protocol _2024_Mar.json`).
"""

from pathlib import Path
import base64
import html
import json
import re
import struct

HERE = Path(__file__).resolve().parent
LUCIDCHART_JSON = HERE / 'Protocol _2024_Mar.json'
BASE_FIGURE = HERE / 'xQTL_Protocol_Feb_2024.png'
OUT = HERE / 'xQTL_Protocol_Implemented_2026.svg'

ADD_ON_COMPONENTS = [
    ('SuSiE-enloc enrichment', 'xQTL-GWAS enrichment priors before pairwise colocalization'),
    ('cTWAS workflow', 'TWAS follow-up with variant selection for candidate genes'),
    ('ColocBoost', 'Additional GWAS/xQTL colocalization method notebook'),
    ('EOO annotation enrichment', 'Excess-of-overlap validation with chromosome jackknife'),
    ('Pathway / GSEA enrichment', 'Gene-set interpretation of prioritized genes'),
    ('GREGOR and sLDSC', 'Variant and partitioned-heritability enrichment workflows'),
    ('xQTL modifier score', 'EMS training and prediction for variant-prioritization support'),
]

EXPECTED_BASE_LABELS = [
    'Existing GWAS & xQTL Summary Statistics',
    'Context Specfic xQTL',
    'Genome-wide xQTL Association',
    'Genome-wide xQTL Fine-mapping',
    'Colocalization[4]',
    'TWAS[8]',
    'Causal TWAS[8]',
    'Mendelian Randomization[4,9]',
    'Multi-Context xQTL',
    "Alzheimer's Disease",
]

COLORS = {
    'panel': '#fff7ed',
    'panel_stroke': '#ea580c',
    'box': '#ffedd5',
    'box_alt': '#e0f2fe',
    'stroke': '#344054',
    'line': '#c2410c',
    'title': '#111827',
    'muted': '#475467',
}


def encoded_png(path):
    return base64.b64encode(path.read_bytes()).decode('ascii')


def png_dimensions(path):
    header = path.read_bytes()[:24]
    if header[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError(f'Expected a PNG file: {path}')
    return struct.unpack('>II', header[16:24])


def lucidchart_labels():
    data = json.loads(LUCIDCHART_JSON.read_text())
    labels = set()
    for page in data.get('pages', []):
        for shape in page.get('items', {}).get('shapes', []):
            text = ' '.join(area.get('text', '') for area in shape.get('textAreas', []))
            text = re.sub(r'\s+', ' ', text).strip()
            if text:
                labels.add(text)
    return labels


def assert_base_figure_is_lucidchart_export():
    labels = lucidchart_labels()
    missing = [label for label in EXPECTED_BASE_LABELS if label not in labels]
    if missing:
        raise ValueError('Expected Lucidchart labels are missing: ' + ', '.join(missing))
    print(f'Audited {len(labels)} labels from {LUCIDCHART_JSON.name}')


def wrap(text, max_chars):
    words = text.split()
    lines, current = [], ''
    for word in words:
        candidate = f'{current} {word}'.strip()
        if len(candidate) <= max_chars:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def text_block(x, y, lines, size, weight='600', fill=None, anchor='middle', line_height=1.2):
    fill = fill or COLORS['title']
    return '\n'.join(
        f'<text x="{x}" y="{y + idx * size * line_height:.1f}" text-anchor="{anchor}" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{weight}" '
        f'fill="{fill}">{html.escape(line)}</text>'
        for idx, line in enumerate(lines)
    )


def box(x, y, w, h, title, subtitle, index):
    fill = COLORS['box_alt'] if index % 2 else COLORS['box']
    title_lines = wrap(title, 24)
    subtitle_lines = wrap(subtitle, 38)
    start_y = y + 74
    parts = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="28" fill="{fill}" stroke="{COLORS["stroke"]}" stroke-width="6"/>']
    parts.append(text_block(x + w / 2, start_y, title_lines, size=64, weight='800'))
    parts.append(text_block(x + w / 2, start_y + 86 * len(title_lines), subtitle_lines, size=44, weight='500', fill=COLORS['muted']))
    return '\n'.join(parts)


def arrow(x1, y1, x2, y2):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{COLORS["line"]}" stroke-width="10" marker-end="url(#arrow)" stroke-dasharray="28 18"/>'


def generate_svg():
    assert_base_figure_is_lucidchart_export()
    base_w, base_h = png_dimensions(BASE_FIGURE)
    add_w = 1680
    gap = 90
    svg_w = base_w + gap + add_w + 90
    svg_h = base_h
    panel_x = base_w + gap
    panel_y = 520
    panel_w = add_w
    panel_h = base_h - 1040
    base_png = encoded_png(BASE_FIGURE)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
        '<defs><marker id="arrow" markerWidth="18" markerHeight="18" refX="15" refY="5" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,10 L16,5 z" fill="#c2410c"/></marker></defs>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<image x="0" y="0" width="{base_w}" height="{base_h}" href="data:image/png;base64,{base_png}"/>',
        f'<rect x="{panel_x}" y="{panel_y}" width="{panel_w}" height="{panel_h}" rx="44" fill="{COLORS["panel"]}" stroke="{COLORS["panel_stroke"]}" stroke-width="10" stroke-dasharray="42 28"/>',
        text_block(panel_x + panel_w / 2, panel_y + 130, ['Repository components to add', 'to the Lucidchart figure'], size=78, weight='900'),
        text_block(panel_x + panel_w / 2, panel_y + 330, ['Base figure is embedded in this SVG;', 'boxes below are missing or newly explicit modules.'], size=44, weight='500', fill=COLORS['muted']),
    ]

    box_x = panel_x + 115
    box_w = panel_w - 230
    box_h = 330
    start_y = panel_y + 560
    step = 420
    for idx, (title, subtitle) in enumerate(ADD_ON_COMPONENTS):
        y = start_y + idx * step
        svg.append(box(box_x, y, box_w, box_h, title, subtitle, idx))
        if idx:
            svg.append(arrow(panel_x + panel_w / 2, y - 92, panel_x + panel_w / 2, y - 8))

    # Dashed callouts from existing integration / validation region into the add-on panel.
    svg.append(arrow(base_w - 260, 2500, panel_x, start_y + 120))
    svg.append(arrow(base_w - 260, 3720, panel_x, start_y + 4 * step + 120))
    svg.append(text_block(panel_x + panel_w / 2, svg_h - 130, [f'Source audit: {LUCIDCHART_JSON.name}; original figure preserved and extended.'], size=42, weight='500', fill=COLORS['muted']))
    svg.append('</svg>')
    OUT.write_text('\n'.join(svg))
    print(OUT)


generate_svg()
