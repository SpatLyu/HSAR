# SAR model estimation

The `sar()` function implements a standard spatial econometrics model
(SAR) or a spatially lagged dependent variable model using the Markov
chain Monte Carlo (McMC) simulation approach.

## Usage

``` r
sar(
  formula,
  data = NULL,
  W,
  Durbin = FALSE,
  burnin = 5000,
  Nsim = 10000,
  thinning = 1,
  parameters.start = NULL
)
```

## Arguments

- formula:

  A symbolic description of the model to fit. A formula for the
  covariate part of the model using the syntax of the
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) function fitting
  standard linear regression models. Neither the response variable nor
  the explanatory variables are allowed to contain NA values.

- data:

  A `data.frame` containing variables used in the formula object.

- W:

  The N by N spatial weights matrix or neighbourhood matrix where N is
  the number of spatial units. The formulation of W could be based on
  geographical distances separating units or based on geographical
  contiguity. To ensure the maximum value of the spatial autoregressive
  parameter \\\rho\\ less than 1, W is usually row-normalised before
  implementing the SAR model. As in most cases, spatial weights matrix
  is very sparse, therefore W here should be converted to a sparse
  matrix before imported into the `sar()` function to save computational
  burden and reduce computing time. More specifically, W should be a
  column-oriented numeric sparse matrices of a `dgCMatrix` class defined
  in the `Matrix` package. The converion between a dense numeric matrix
  and a sparse numeric matrix is made quite convenient through the
  `Matrix` library.

- Durbin:

  `logical`. Estimate Durbin model (i.e. include spatial lags of `X` as
  predictors)? Default `FALSE`.

- burnin:

  The number of McMC samples to discard as the burnin period.

- Nsim:

  The total number of McMC samples to generate.

- thinning:

  MCMC thinning factor.

- parameters.start:

  A list with names "rho", "sigma2e", and "beta" corresponding to
  initial values for the model parameters \\\rho, \sigma^2_e\\ and the
  regression coefficients respectively.

## Value

A `list`.

- cbetas:

  A matrix with the MCMC samples of the draws for the coefficients.

- Mbetas:

  A vector of estimated mean values of regression coefficients.

- SDbetas:

  The standard deviations of estimated regression coefficients.

- Mrho:

  The estimated mean of the lower-level spatial autoregressive parameter
  \\\rho\\.

- SDrho:

  The standard deviation of the estimated lower-level spatial
  autoregressive parameter.

- Msigma2e:

  The estimated mean of the lower-level variance parameter
  \\\sigma^{2}\_{e} \\.

- SDsigma2e:

  The standard deviation of the estimated lower-level variance parameter
  \\\sigma^{2}\_{e} \\.

- DIC:

  The deviance information criterion (DIC) of the fitted model.

- pd:

  The effective number of parameters of the fitted model.

- Log_Likelihood:

  The log-likelihood of the fitted model.

- R_Squared:

  A pseudo R square model fit indicator.

- impact_direct:

  Summaries of the direct impact of a covariate effect on the outcome
  variable.

- impact_idirect:

  Summaries of the indirect impact of a covariate effect on the outcome
  variable.

- impact_total:

  Summaries of the total impact of a covariate effect on the outcome
  variable.

## References

Anselin, L. (1988). *Spatial Econometrics: Methods and Models*.
Dordrecht: Kluwer Academic Publishers.

LeSage, J. P., and R. K. Pace. (2009). *Introduction to Spatial
Econometrics*. Boca Raton, FL: CRC Press/Taylor & Francis

## Examples

``` r
data(landprice)
head(landprice)
#>   obs lnprice   lnarea  lndcbd dsubway   dpark    dele   popden crimerate
#> 1 190 7.16382 11.58780 9.93534 7.14334 6.78243 6.67827 0.548966   10.7511
#> 2 992 7.45757  9.20029 9.84785 7.78904 6.95662 7.05138 0.548966   10.7511
#> 3 189 5.57430 10.27820 9.94866 6.83023 7.06579 6.81916 0.548966   10.7511
#> 4 993 7.12569  7.97788 9.84388 7.81991 7.00792 7.11267 0.548966   10.7511
#> 5 969 6.81564  5.81928 9.91940 7.64640 6.88254 4.10025 0.548966   10.7511
#> 6 994 7.48522  7.78634 9.84203 7.83398 7.03089 7.13981 0.548966   10.7511
#>   district.id year
#> 1           3    1
#> 2           3    0
#> 3           3    1
#> 4           3    0
#> 5           3    0
#> 6           3    0
data(land)

# extract the land parcel level spatial weights matrix
library(spdep)
library(Matrix)
nb.25 <- spdep::dnearneigh(land,0,2500)
#> Warning: neighbour object has 4 sub-graphs
# to a weights matrix
dist.25 <- spdep::nbdists(nb.25,land)
dist.25 <- lapply(dist.25,function(x) exp(-0.5 * (x / 2500)^2))
mat.25 <- spdep::nb2mat(nb.25,glist=dist.25,style="W")
W <- as(mat.25,"dgCMatrix")

## run the sar() function
res.formula <- lnprice ~ lnarea + lndcbd + dsubway + dpark + dele +
                popden + crimerate + as.factor(year)
betas= coef(lm(formula=res.formula,data=landprice))
pars=list(rho = 0.5, sigma2e = 2.0, betas = betas)
# \donttest{
res <- sar(res.formula,data=landprice,W=W,
           burnin=500, Nsim=1000, thinning=1,
           parameters.start=pars)
#> Warning: style is M (missing); style should be set to a valid value
#> Warning: style is M (missing); style should be set to a valid value
#> Warning: neighbour object has 4 sub-graphs
#> Warning: neighbour object has 4 sub-graphs
summary(res)
#> 
#> Call:
#> sar(formula = res.formula, data = landprice, W = W, burnin = 500, 
#>     Nsim = 1000, thinning = 1, parameters.start = pars)
#> Type:  sar  
#> 
#> Coefficients:
#>                          Mean          SD
#> (Intercept)       6.692659271 0.821582708
#> lnarea           -0.005703754 0.017245283
#> lndcbd           -0.093812794 0.048795925
#> dsubway          -0.151066147 0.035977010
#> dpark            -0.157082432 0.045501847
#> dele             -0.029478755 0.031136321
#> popden            0.021820387 0.010123086
#> crimerate         0.006581084 0.004185754
#> as.factor(year)1 -0.188791418 0.055353624
#> as.factor(year)2 -0.028732593 0.117591116
#> as.factor(year)3 -0.141198327 0.106523962
#> as.factor(year)4  0.568283074 0.116757966
#> as.factor(year)5  0.398618141 0.127220897
#> as.factor(year)6  2.183383251 0.218838610
#> 
#>  Spatial Coefficients:
#>      rho 
#> 0.544068 
#> 
#>  Diagnostics 
#> Deviance information criterion (DIC): 4704.887 
#> Effective number of parameters (pd): 956.0725 
#> Log likelihood: -1396.371 
#> Pseudo R squared: 0.3487775 
#> 
#>  Impacts:
#>                        direct     indirect       total
#> lnarea           -0.005808262 -0.006377362 -0.01218562
#> lndcbd           -0.095531697 -0.104891996 -0.20042369
#> dsubway          -0.153834086 -0.168906916 -0.32274100
#> dpark            -0.159960606 -0.175633718 -0.33559432
#> dele             -0.030018885 -0.032960168 -0.06297905
#> popden            0.022220195  0.024397354  0.04661755
#> crimerate         0.006701667  0.007358304  0.01405997
#> as.factor(year)1 -0.192250585 -0.211087505 -0.40333809
#> as.factor(year)2 -0.029259052 -0.032125885 -0.06138494
#> as.factor(year)3 -0.143785461 -0.157873715 -0.30165918
#> as.factor(year)4  0.578695551  0.635396766  1.21409232
#> as.factor(year)5  0.405921899  0.445694565  0.85161646
#> as.factor(year)6  2.223388716  2.441238743  4.66462746
#> 
#>  Quantiles:
#>                             5%          25%          50%          75%
#> (Intercept)       5.3931321636  6.102344739  6.667925879  7.261131168
#> lnarea           -0.0328579787 -0.017805259 -0.005443984  0.005865107
#> lndcbd           -0.1719583291 -0.126473435 -0.094248991 -0.061389837
#> dsubway          -0.2086087870 -0.174525097 -0.151688345 -0.127059403
#> dpark            -0.2366361542 -0.187072950 -0.157027536 -0.126171278
#> dele             -0.0799836946 -0.051391759 -0.029007668 -0.008189104
#> popden            0.0054129196  0.014553396  0.021759805  0.028808175
#> crimerate        -0.0004581391  0.003773817  0.006547580  0.009282714
#> as.factor(year)1 -0.2756505015 -0.227269916 -0.190331369 -0.150788147
#> as.factor(year)2 -0.2218306150 -0.106403948 -0.031509779  0.056263472
#> as.factor(year)3 -0.3146393375 -0.212914997 -0.142246504 -0.070137649
#> as.factor(year)4  0.3816845823  0.487045899  0.577141855  0.645571177
#> as.factor(year)5  0.1922468443  0.311307645  0.402879081  0.481189056
#> as.factor(year)6  1.8106568356  2.041326641  2.191568289  2.327773608
#>                          95%
#> (Intercept)       8.03383299
#> lnarea            0.02201095
#> lndcbd           -0.01651091
#> dsubway          -0.09292392
#> dpark            -0.08466489
#> dele              0.02294826
#> popden            0.03870791
#> crimerate         0.01327714
#> as.factor(year)1 -0.09626197
#> as.factor(year)2  0.15447486
#> as.factor(year)3  0.04359984
#> as.factor(year)4  0.76317569
#> as.factor(year)5  0.60496796
#> as.factor(year)6  2.55235508
# }
```
