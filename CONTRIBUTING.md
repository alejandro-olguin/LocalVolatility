# Contributing

Thanks for your interest in improving LocalVolatility!

- Use feature branches and pull requests.
- Add tests for new functionality and keep coverage high.
- Run R CMD check locally (or rely on CI) and fix `lintr` issues before submitting.
- Follow the existing coding style (Rcpp C++17, roxygen2 docs colocated with exports).

## Dev quickstart

```r
Rcpp::compileAttributes()
roxygen2::roxygenise()
devtools::load_all()
devtools::test()
```

Please open issues at https://github.com/alejandro-olguin/LocalVolatility/issues.
