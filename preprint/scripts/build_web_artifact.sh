#!/usr/bin/env bash
# Build the canonical PDF rendering of the preprint and vendor
# it into ldsc-web/assets/preprint.pdf, where the Preprint tab
# iframes it directly.
#
# Output (committed; CI doesn't need typst):
#   ldsc-web/assets/preprint.pdf — high-fidelity 2 MB PDF, rendered
#                                  by the browser's built-in PDF
#                                  viewer (PDFium / pdf.js / PDFKit).
#                                  No JS shim, no SVG hackery —
#                                  selection / Ctrl-F / zoom all
#                                  come from the browser's reader.
#
# Re-run whenever preprint/main.typ, preprint/template.typ, or any
# preprint/figures/*.png changes, then commit the regenerated PDF.
#
# Producer requirement (one-time install):
#   brew install typst    # or `cargo install --locked typst-cli`

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

if ! command -v typst >/dev/null 2>&1; then
    echo "error: typst not on PATH" >&2
    echo "install: brew install typst   # or: cargo install --locked typst-cli" >&2
    exit 1
fi

echo "==> compiling preprint/main.typ → PDF"
typst compile preprint/main.typ preprint/main.pdf

echo "==> vendoring → ldsc-web/assets/preprint.pdf"
mkdir -p ldsc-web/assets
cp preprint/main.pdf ldsc-web/assets/preprint.pdf

echo
echo "vendored artifact:"
ls -lh ldsc-web/assets/preprint.pdf
