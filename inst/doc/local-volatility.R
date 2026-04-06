## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----eval=FALSE---------------------------------------------------------------
#  library(LocalVolatility)
#  # Quick constant-vol European 2D price
#  price <- european_option_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
#                              0.2, 0.15, 0.3, "call",
#                              10, 300, 1, 60, 40, 40, 50, 3)
#  price
