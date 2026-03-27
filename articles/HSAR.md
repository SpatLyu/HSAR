# Hierarchical Spatial Simultaneous Autoregressive Model (HSAR)

## An application of HSAR for asking prices in the municipality of Athens

An application of
[`hsar()`](https://spatlyu.github.io/HSAR/reference/hsar.md), based on
rel data, will be illustrated. The design of the weight matrices needed
and the random effect design matrix will be explained.

#### Libraries

We start by loading the libraries that will be used.

``` r
library(sf)
library(spdep)
library(tidyverse)
## ── Attaching core tidyverse packages ─────────────────────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.0.2     
## ── Conflicts ───────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ tidyr::expand() masks Matrix::expand()
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ✖ tidyr::pack()   masks Matrix::pack()
## ✖ tidyr::unpack() masks Matrix::unpack()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(HSAR)
```

#### Reading the datasets

At the higher level, we have the seven departments of the municipality
of Athens and at the lower level we have the point data of the
properties.

``` r
data(depmunic)
data(properties)
plot(st_geometry(depmunic),col = sf.colors(12, categorical = TRUE), border = 'grey')
plot(st_geometry(properties),add=TRUE,col="red",pch=16,cex=0.6)
```

![](../reference/figures/hsar/p1-1.png)

The characteristics that come with the areal data are the id of the
department, the number of airbnb properties, the number of museums, the
population, the number of citizens with origin a non european union
country, the area of the green space (m^2) and the area of the polygon
(km^2).

``` r
names(depmunic)
## [1] "num_dep"    "airbnb"     "museums"    "population" "pop_rest"   "greensp"    "area"      
## [8] "geometry"
depmunic$pop_rest
## [1]  8202  5009  2735  4167  5099 16531  8017
```

The characteristics of the properties are the size (m^2), the asking
price (euros), the price per square meter, the age (years) and the
shortest distance to metro/train station (m).

``` r
names(properties)
## [1] "id"         "size"       "price"      "prpsqm"     "age"        "dist_metro" "geometry"
hist(properties$age, xlab = "Age", main="Age of the properties")
```

![](../reference/figures/hsar/p2-1.png)

Now we are going to create two more variables at the higher,
municipality department, level. The first one is the population density
per 10k citizens, and the second one is the percentage of non EU
citizens.

``` r
depmunic$popdens <- depmunic$population/ (10000*depmunic$area)
depmunic$foreigners <- 100 * depmunic$pop_rest/ depmunic$population
```

The next step is to create the model data that are going to use in the
hsar model. For that, we need for each property (lower data), the data
from the relevant department(higher level).

``` r
properties_in_dd <- st_join(properties, depmunic, join = st_within)
```

So now, we know each property, in which department resides and the
coresponding data for that polygon. We also need that data in sorting
order.

``` r
model.data <- properties_in_dd[order(properties_in_dd$num_dep),]
```

#### Create matrices used in the hsar function

In order to run the model we need to create the effect design matrix
(Delta), the weight matrix for the high-level - polygon data (M), and
the weight matrix for the lower level - point data (W).

In order to define the random effect matrix, we start with estimating
the number of properties in each municipality department

``` r
properties_count <- count(as_tibble(model.data), num_dep)
MM <- as.data.frame(properties_count)
```

and by geting the total number of municipality departments (7), we
define a vector with the number of municipality department that each
property belongs

``` r
Utotal <- dim(MM)[1]
Unum <- MM[,2]
Uid <- rep(c(1:Utotal),Unum)
```

We then define the random effect matrix (Delta) wich has a dimension of
1000x7

``` r
n <- nrow(properties)
Delta <- matrix(0,nrow=n,ncol=Utotal)
for(i in 1:Utotal) {
  Delta[Uid==i,i] <- 1
}

Delta <- as(Delta,"dgCMatrix")
```

Now we estimate the spatial weight matrix at the higher level which in
our case is the municipality departments (polygons). So we start with
poly2nb which constructs the neighbours list for polygons and then with
nb2mat we generate the weight matrix for the neighbours list previously
created. Then we transform the weight matrix in a sparse matrix format.

``` r
nb.list <- poly2nb(depmunic)
mat.list <- nb2mat(nb.list,style="W")
M <- as(mat.list,"dgCMatrix")
```

to have a closer look at M , we can visualize it

``` r
plot(st_geometry(depmunic),border = 'grey')
plot(st_centroid(depmunic), add = TRUE)
## Warning: st_centroid assumes attributes are constant over geometries
## Warning in plot.sf(st_centroid(depmunic), add = TRUE): ignoring all but the first attribute
plot(nb.list, st_centroid(depmunic), add = TRUE)
## Warning: st_centroid assumes attributes are constant over geometries
```

![](../reference/figures/hsar/p3-1.png)

Similarly, we create the spatial weight matrix at the lower level of
properties (point data). So we create the neighbour list at a distance
of 1300 meters

``` r
nb.1300 <- dnearneigh(properties,0,1300)
```

and the weights matrix W as follows

``` r
mat.1300 <- nb2mat(nb.1300,style="W")
W <- as(mat.1300,"dgCMatrix")
```

For the W matrix, we can check the neighbours statistics

``` r
nb.1300
## Neighbour list object:
## Number of regions: 1000 
## Number of nonzero links: 170254 
## Percentage nonzero weights: 17.0254 
## Average number of links: 170.254
```

#### Run the models

So, having ready the matrices Delta, M and W, we wun the
[`hsar()`](https://spatlyu.github.io/HSAR/reference/hsar.md) function

``` r
res.formula <- prpsqm ~ size + age + greensp + population + museums + airbnb
res <- hsar(res.formula,data=model.data,W=W,M=M,Delta=Delta,
            burnin=500, Nsim=1000)
## Warning in spdep::mat2listw(W): style is M (missing); style should be set to a valid value
## Warning in sn2listw(df, style = style, zero.policy = zero.policy, from_mat2listw = TRUE):
## style is M (missing); style should be set to a valid value
## Warning in spdep::mat2listw(W): style is M (missing); style should be set to a valid value
## Warning in sn2listw(df, style = style, zero.policy = zero.policy, from_mat2listw = TRUE):
## style is M (missing); style should be set to a valid value
summary(res)
## 
## Call:
## hsar(formula = res.formula, data = model.data, W = W, M = M, 
##     Delta = Delta, burnin = 500, Nsim = 1000)
## Type:  hsar  
## 
## Coefficients:
##                      Mean           SD
## (Intercept)  1.880468e+03 9.835447e+00
## size         4.298802e+00 5.112041e-01
## age         -1.995687e+01 1.304370e+00
## greensp      8.404794e-04 8.738884e-04
## population  -9.940391e-03 2.273935e-03
## museums     -4.515772e+01 1.008513e+01
## airbnb       6.022747e-01 2.496115e-01
## 
##  Spatial Coefficients:
##           rho   lambda
## [1,] 0.196536 0.018432
## 
##  Diagnostics 
## Deviance information criterion (DIC): 28193.01 
## Effective number of parameters (pd): -1.66553 
## Log likelihood: -14098.17 
## Pseudo R squared: 0.3601049 
## 
##  Impacts:
##                    direct      indirect         total
## (Intercept)  1.881082e+03  4.592335e+02  2.340316e+03
## size         4.300207e+00  1.049821e+00  5.350028e+00
## age         -1.996339e+01 -4.873715e+00 -2.483711e+01
## greensp      8.407540e-04  2.052555e-04  1.046009e-03
## population  -9.943639e-03 -2.427567e-03 -1.237121e-02
## museums     -4.517248e+01 -1.102808e+01 -5.620055e+01
## airbnb       6.024715e-01  1.470830e-01  7.495545e-01
## 
##  Quantiles:
##                        5%           25%           50%           75%           95%
## (Intercept)  1.864369e+03  1.874059e+03  1.880361e+03  1.887328e+03  1.896080e+03
## size         3.436255e+00  3.956910e+00  4.291712e+00  4.628642e+00  5.155081e+00
## age         -2.215863e+01 -2.075471e+01 -1.991798e+01 -1.911688e+01 -1.792525e+01
## greensp     -3.546526e-04  2.905312e-04  7.376810e-04  1.261277e-03  2.657171e-03
## population  -1.382163e-02 -1.159143e-02 -9.920309e-03 -8.353503e-03 -6.404121e-03
## museums     -6.148110e+01 -5.224350e+01 -4.495013e+01 -3.850940e+01 -2.847943e+01
## airbnb       1.776907e-01  4.578970e-01  6.114716e-01  7.457536e-01  1.010643e+00
```

and the two simpler models defined for rho = 0 and lambda=0. So,
firstly, assuming rho = 0 (no interaction effects at the lower level) we
get

``` r
res_1 <- hsar(res.formula,data=model.data,W=NULL,M=M,Delta=Delta,burnin=500, Nsim=1000)
## Warning in spdep::mat2listw(W): style is M (missing); style should be set to a valid value
## Warning in sn2listw(df, style = style, zero.policy = zero.policy, from_mat2listw = TRUE):
## style is M (missing); style should be set to a valid value
summary(res_1)
## 
## Call:
## hsar(formula = res.formula, data = model.data, W = NULL, M = M, 
##     Delta = Delta, burnin = 500, Nsim = 1000)
## Type:  hsar with rho = 0  
## 
## Coefficients:
##                      Mean           SD
## (Intercept)  1.880592e+03 1.010450e+01
## size         4.329009e+00 4.214029e-01
## age         -2.004845e+01 1.258780e+00
## greensp      6.812935e-04 6.843458e-04
## population  -6.877688e-03 1.100394e-03
## museums     -4.585817e+01 9.479326e+00
## airbnb       6.288388e-01 2.222552e-01
## 
##  Spatial Coefficients:
##    lambda 
## -0.130888 
## 
##  Diagnostics 
## Deviance information criterion (DIC): 28196.24 
## Effective number of parameters (pd): -1.988165 
## Log likelihood: -14100.11 
## Pseudo R squared: 0.3587056 
## 
##  Quantiles:
##                        5%           25%           50%           75%           95%
## (Intercept)  1.863327e+03  1.874280e+03  1.880446e+03  1.887451e+03  1.897526e+03
## size         3.648867e+00  4.026734e+00  4.331703e+00  4.647616e+00  5.040309e+00
## age         -2.217344e+01 -2.090501e+01 -2.002926e+01 -1.916105e+01 -1.800675e+01
## greensp     -5.464395e-04  3.468840e-04  7.498473e-04  1.081365e-03  1.660609e-03
## population  -8.718853e-03 -7.542340e-03 -6.897453e-03 -6.200430e-03 -4.988964e-03
## museums     -6.078241e+01 -5.205170e+01 -4.542701e+01 -3.995038e+01 -3.145765e+01
## airbnb       2.844778e-01  4.953275e-01  6.164306e-01  7.496283e-01  9.872352e-01
```

and secondly, given lambda = 0 (no interaction at the higher level) we
get

``` r
res_2 <- hsar(res.formula,data=model.data,W=W,M=NULL,Delta=Delta,burnin=500, Nsim=1000)
## Warning in spdep::mat2listw(W): style is M (missing); style should be set to a valid value
## Warning in sn2listw(df, style = style, zero.policy = zero.policy, from_mat2listw = TRUE):
## style is M (missing); style should be set to a valid value
summary(res_2)
## 
## Call:
## hsar(formula = res.formula, data = model.data, W = W, M = NULL, 
##     Delta = Delta, burnin = 500, Nsim = 1000)
## Type:  hsar with lambda = 0  
## 
## Coefficients:
##                      Mean           SD
## (Intercept)  1.880293e+03 9.970170e+00
## size         4.271588e+00 4.581010e-01
## age         -1.995568e+01 1.295576e+00
## greensp      9.672711e-04 6.581401e-04
## population  -9.434666e-03 2.083913e-03
## museums     -4.545286e+01 1.035321e+01
## airbnb       5.415177e-01 1.964108e-01
## 
##  Spatial Coefficients:
##     rho 
## 0.19122 
## 
##  Diagnostics 
## Deviance information criterion (DIC): 28196.99 
## Effective number of parameters (pd): -1.752399 
## Log likelihood: -14100.25 
## Pseudo R squared: 0.3597633 
## 
##  Quantiles:
##                        5%           25%           50%           75%           95%
## (Intercept)  1.864332e+03  1.873364e+03  1.880298e+03  1.886796e+03  1.896416e+03
## size         3.550689e+00  3.944261e+00  4.254099e+00  4.559823e+00  5.073155e+00
## age         -2.192657e+01 -2.079846e+01 -1.995871e+01 -1.909676e+01 -1.780104e+01
## greensp      3.886899e-05  5.256176e-04  8.835621e-04  1.378462e-03  2.197482e-03
## population  -1.282486e-02 -1.074985e-02 -9.589800e-03 -8.107401e-03 -6.107883e-03
## museums     -6.271483e+01 -5.293104e+01 -4.448211e+01 -3.801137e+01 -2.948482e+01
## airbnb       2.208243e-01  4.100049e-01  5.411334e-01  6.705060e-01  8.635070e-01
```
