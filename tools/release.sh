#!/usr/bin/env bash
set -euo pipefail

ver=${1:-2.0.0}

R -q -e "Rcpp::compileAttributes(); roxygen2::roxygenise(); devtools::check(document = FALSE, manual = TRUE)"

git add -A
if ! git diff --cached --quiet; then
  git commit -m "Release ${ver}"
fi

git tag -a "v${ver}" -m "LocalVolatility ${ver}" || true
git push && git push --tags
