
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HSAR

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/HSAR)](https://CRAN.R-project.org/package=HSAR)
[![DownloadsAll](https://badgen.net/cran/dt/HSAR?color=orange)](https://CRAN.R-project.org/package=HSAR)
[![DownloadsMonth](https://cranlogs.r-pkg.org/badges/HSAR)](https://CRAN.R-project.org/package=HSAR)
[![R-universe](https://spatlyu.r-universe.dev/badges/HSAR?color=cyan)](https://spatlyu.r-universe.dev/HSAR)
[![R-CMD-check](https://github.com/SpatLyu/HSAR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SpatLyu/HSAR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Hierarchical Spatial Simultaneous Autoregressive Model (HSAR)**

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=HSAR) with:

``` r
install.packages("HSAR", dep = TRUE)
```

- Install development binary version from
  [R-universe](https://spatlyu.r-universe.dev/HSAR) with:

``` r
install.packages('HSAR',
                 repos = c("https://spatlyu.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install development source version from
  [GitHub](https://github.com/spatlyu/HSAR) with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("SpatLyu/HSAR",
                         build_vignettes = TRUE,
                         dep = TRUE)
```
