pkgname <- "LocalVolatility"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('LocalVolatility')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("american_option_heston")
### * american_option_heston

flush(stderr()); flush(stdout())

### Name: american_option_heston
### Title: American Option (Heston stochastic volatility, penalty method)
### Aliases: american_option_heston

### ** Examples

american_option_heston(100, 100, 1, 0.05, 0, 2, 0.04, 0.5, -0.7, 0.04,
                       "put", 20, 300, 0.001, 1.0, 80, 40, 100, 3,
                       1e4, 1e-8)




cleanEx()
nameEx("european_option_cf_2d")
### * european_option_cf_2d

flush(stderr()); flush(stdout())

### Name: european_option_cf_2d
### Title: European Option 2D (closed-form, ADR)
### Aliases: european_option_cf_2d

### ** Examples

european_option_cf_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01, 0.2, 0.1, 0.3, "call", 1)




cleanEx()
nameEx("european_option_heston")
### * european_option_heston

flush(stderr()); flush(stdout())

### Name: european_option_heston
### Title: European Option (Heston stochastic volatility)
### Aliases: european_option_heston

### ** Examples

european_option_heston(100, 100, 1, 0.05, 0, 2, 0.04, 0.5, -0.7, 0.04,
                       "call", 20, 300, 0.001, 1.0, 80, 40, 100, 3)




cleanEx()
nameEx("heston_cf")
### * heston_cf

flush(stderr()); flush(stdout())

### Name: heston_cf
### Title: Heston Model (closed-form, characteristic function)
### Aliases: heston_cf

### ** Examples

heston_cf(100, 100, 1, 0.05, 0, 2, 0.04, 0.5, -0.7, 0.04, "call")




cleanEx()
nameEx("mc_american_heston_4d")
### * mc_american_heston_4d

flush(stderr()); flush(stdout())

### Name: mc_american_heston_4d
### Title: American Option 4D Monte Carlo (Double-Heston Quanto,
###   Longstaff-Schwartz)
### Aliases: mc_american_heston_4d

### ** Examples

mc_american_heston_4d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
  2, 0.04, 0.5, 0.04, 1.5, 0.02, 0.3, 0.02,
  0.3, -0.7, -0.5, 0, 0, 0, "put", 50000, 50, 42, 4)




cleanEx()
nameEx("mc_european_heston_4d")
### * mc_european_heston_4d

flush(stderr()); flush(stdout())

### Name: mc_european_heston_4d
### Title: European Option 4D Monte Carlo (Double-Heston Quanto)
### Aliases: mc_european_heston_4d

### ** Examples

mc_european_heston_4d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
  2, 0.04, 0.5, 0.04, 1.5, 0.02, 0.3, 0.02,
  0.3, -0.7, -0.5, 0, 0, 0, "call", 100000, 100, 42)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
