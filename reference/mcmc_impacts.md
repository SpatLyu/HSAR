# Estimate impact per MCMC draw of coefficients

Estimate impact per MCMC draw of coefficients

## Usage

``` r
mcmc_impacts(object, p = 5)
```

## Arguments

- object:

  an `mcmc_hsar` or `mcmc_hsar_lambda_0` or `mcmc_sar` object

- p:

  `numeric`. Maximum weight matrix power. Default `5`.

## Value

a list of three matrices named `direct`, `indirect`, and `total` giving
impact estimates for each MCMC draw
