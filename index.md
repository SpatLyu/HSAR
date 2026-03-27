# HSAR

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
