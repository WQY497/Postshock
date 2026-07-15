#!/bin/zsh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
R_SCRIPT="/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"
PANDOC_DIR="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"

cd "$ROOT_DIR"

RSTUDIO_PANDOC="$PANDOC_DIR" \
"$R_SCRIPT" -e 'rmarkdown::render("inst/technical_report/Technical_report.Rmd", output_format = "pdf_document")'
