#!/usr/bin/env bash
# Build the standalone HTML rendering of the preprint that the
# ldsc-web "Preprint" tab iframes.
#
# Outputs (committed; CI doesn't need typst-ts-cli):
#   ldsc-web/assets/preprint.html — self-contained HTML with the
#                                    full preprint as a single
#                                    inline SVG (selectable text
#                                    via .tsel overlay, all CSS
#                                    inlined). ~8 MB raw, ~2.3 MB
#                                    gzipped on the wire.
#   ldsc-web/assets/preprint.pdf  — high-fidelity PDF download
#                                    fallback (2 MB).
#
# Re-run whenever preprint/main.typ, preprint/template.typ, or any
# preprint/figures/*.png changes, then commit the regenerated
# outputs.
#
# One-time install of the producer:
#   cargo install --locked --git https://github.com/Myriad-Dreamin/typst.ts typst-ts-cli
# (We use Myriad-Dreamin's fork of typst, which adds the `svg_html`
# format — upstream typst doesn't have it. Version 0.7.0 confirmed
# working.)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

if ! command -v typst-ts-cli >/dev/null 2>&1; then
    echo "error: typst-ts-cli not on PATH" >&2
    echo "install: cargo install --locked --git https://github.com/Myriad-Dreamin/typst.ts typst-ts-cli" >&2
    exit 1
fi

mkdir -p ldsc-web/assets

# Compile into a temp dir, then post-process into the final
# location. The tmp step is so we don't leave a half-baked
# preprint.html sitting in assets/ if either step fails.
TMP="$(mktemp -d -t ldsc-preprint-build.XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

echo "==> compiling preprint → svg_html"
# `--workspace preprint` makes typst-ts-cli treat the preprint
# directory as the project root (so `image("figures/...")` paths
# inside main.typ resolve correctly).
# `--format svg_html` emits a single standalone HTML doc with the
# preprint rendered as one tall inline SVG, plus the typst.ts
# semantic-text overlay (`.tsel` class) for native selection /
# Ctrl-F search.
typst-ts-cli compile \
    --workspace preprint \
    --entry preprint/main.typ \
    --format svg_html \
    --output "$TMP/"

# typst-ts-cli emits with the fixed name `main.artifact.svg.html`.
SRC="$TMP/main.artifact.svg.html"
if [[ ! -f "$SRC" ]]; then
    echo "error: typst-ts-cli did not produce $SRC" >&2
    exit 1
fi

# Post-process: the root <svg> element ships with hard-coded
# `width="612.000" height="16632.000"` attributes (US-letter
# width × 20 pages tall). Strip those so the SVG scales
# responsively to the iframe container width — critical for
# mobile, where 612 px would force horizontal scroll and clip
# content. The viewBox is preserved (controls the aspect ratio)
# and we inject `style="display:block; width:100%; height:auto;"`
# so the SVG fills the iframe width and the iframe's own height
# constraint (set in palette.css) controls vertical sizing.
echo "==> post-processing for responsive width"
python3 - <<PY
import re
src = "$SRC"
dst = "ldsc-web/assets/preprint.html"
content = open(src, encoding="utf-8").read()

# Match the first <svg ...> tag (the document root). Use a single
# replacement to (a) drop width/height attrs, (b) inject the
# responsive style. Use a function so we can preserve other attrs
# (viewBox, class, xmlns, etc.) verbatim.
def fix_svg(m):
    attrs = m.group(1)
    attrs = re.sub(r'\s+(?:width|height)="[\d.]+"', '', attrs)
    attrs = re.sub(r'\s+style="[^"]*"', '', attrs)
    return '<svg' + attrs + ' style="display:block; width:100%; height:auto;">'

content, n = re.subn(r'<svg([^>]+?)>', fix_svg, content, count=1)
if n != 1:
    raise SystemExit(f"expected exactly 1 root <svg>, replaced {n}")

open(dst, "w", encoding="utf-8").write(content)
PY

# Vendor the PDF as a download fallback. Trunk's copy-file
# preserves source filenames so we land the copy here with the
# final dist/ name we want.
if [[ -f preprint/main.pdf ]]; then
    cp preprint/main.pdf ldsc-web/assets/preprint.pdf
else
    echo "warning: preprint/main.pdf not found — skipping PDF copy" >&2
    echo "         (rebuild with: typst compile preprint/main.typ)" >&2
fi

echo
echo "vendored artifacts:"
ls -lh ldsc-web/assets/preprint.html ldsc-web/assets/preprint.pdf 2>/dev/null
echo
echo "gzipped size estimate (what GH Pages serves):"
gzip -c ldsc-web/assets/preprint.html | wc -c | awk '{printf "  preprint.html: %.2f MB\n", $1/1024/1024}'
