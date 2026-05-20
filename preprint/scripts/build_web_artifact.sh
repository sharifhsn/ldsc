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

# Post-process two things:
#
# (1) The root <svg> element ships with hard-coded
#     `width="612.000" height="16632.000"` attributes (US-letter
#     width × 20 pages tall). Strip those so the SVG scales
#     responsively to the iframe container width — critical for
#     mobile, where 612 px would force horizontal scroll and clip
#     content. The viewBox is preserved (controls the aspect ratio)
#     and we inject `style="display:block; width:100%; height:auto;"`
#     so the SVG fills the iframe width and the iframe's own height
#     constraint (set in palette.css) controls vertical sizing.
#
# (2) Inject a tiny drag-selection shim before </body>. typst.ts
#     svg_html wraps every text run in its own <foreignObject>
#     (6,946 of them in this preprint), and browsers' native
#     drag-selection refuses to extend a Range across SVG
#     foreignObject boundaries even though caretPositionFromPoint
#     works fine on each individually — so click-and-drag produces
#     an empty selection unless the user double-clicks first.
#     The shim takes over drag-selection: mousedown records a
#     caret, mousemove rebuilds the Range via
#     getSelection().setBaseAndExtent() targeting whatever caret
#     the cursor is currently over. Double/triple-click are NOT
#     intercepted; the browser's word/paragraph-selection logic
#     still handles those natively. ~40 LOC, ~1.6 KB inline.
#     Requires `sandbox="allow-scripts"` on the iframe in
#     preprint_panel.rs.
echo "==> post-processing for responsive width + selection shim"
SRC="$SRC" python3 - <<'PY'
import os
import re
src = os.environ["SRC"]
dst = "ldsc-web/assets/preprint.html"
content = open(src, encoding="utf-8").read()

# (1) Responsive root <svg>.
def fix_svg(m):
    attrs = m.group(1)
    attrs = re.sub(r'\s+(?:width|height)="[\d.]+"', '', attrs)
    attrs = re.sub(r'\s+style="[^"]*"', '', attrs)
    return '<svg' + attrs + ' style="display:block; width:100%; height:auto;">'

content, n = re.subn(r'<svg([^>]+?)>', fix_svg, content, count=1)
if n != 1:
    raise SystemExit(f"expected exactly 1 root <svg>, replaced {n}")

# (2) Drag-selection shim, injected before </body>.
shim = """<script>
(function(){
  // Drag-selection shim for typst.ts svg_html.
  // Each text run is in its own <foreignObject>; browsers refuse
  // to extend a Range across that boundary during mouse drag.
  // We drive selection manually via caretPositionFromPoint.
  function caretAt(x, y){
    if (document.caretPositionFromPoint){
      var c = document.caretPositionFromPoint(x, y);
      return c ? {n: c.offsetNode, o: c.offset} : null;
    }
    if (document.caretRangeFromPoint){
      var r = document.caretRangeFromPoint(x, y);
      return r ? {n: r.startContainer, o: r.startOffset} : null;
    }
    return null;
  }
  var start = null;
  document.addEventListener('mousedown', function(e){
    if (e.button !== 0) return;
    var t = e.target;
    if (t && t.classList && t.classList.contains('pseudo-link')) return;
    var c = caretAt(e.clientX, e.clientY);
    if (!c) { start = null; return; }
    start = c;
    var s = window.getSelection();
    if (s) s.removeAllRanges();
  }, true);
  document.addEventListener('mousemove', function(e){
    if (!start) return;
    if (!(e.buttons & 1)) { start = null; return; }
    var end = caretAt(e.clientX, e.clientY);
    if (!end) return;
    var s = window.getSelection();
    try { s.setBaseAndExtent(start.n, start.o, end.n, end.o); e.preventDefault(); }
    catch (_) { /* invalid range; ignore */ }
  }, true);
  document.addEventListener('mouseup', function(){ start = null; }, true);
})();
</script>"""

if "</body>" in content:
    content = content.replace("</body>", shim + "</body>", 1)
else:
    # Some typst.ts builds emit no explicit </body>; append at end.
    content += shim

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
