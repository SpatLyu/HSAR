#' SAR model estimation
#' @description
#' The `sar()` function implements a standard spatial econometrics model (SAR) or a spatially lagged dependent
#' variable model using the Markov chain Monte Carlo (McMC) simulation approach.
#' @references
#' Anselin, L. (1988). \emph{Spatial Econometrics: Methods and Models}. Dordrecht: Kluwer Academic Publishers.
#'
#' LeSage, J. P., and R. K. Pace. (2009). \emph{Introduction to Spatial Econometrics}. Boca Raton, FL: CRC Press/Taylor & Francis
#'
#' @param formula A symbolic description of the model to fit. A formula for the covariate part of the model
#' using the syntax of the `stats::lm()` function fitting standard linear regression models. Neither the
#' response variable nor the explanatory variables are allowed to contain NA values.
#' @param data A `data.frame` containing variables used in the formula object.
#' @param W The N by N spatial weights matrix or neighbourhood matrix where N is the number of spatial units.
#' The formulation of W could be based on geographical distances separating units or based on geographical contiguity.
#' To ensure the maximum value of the spatial autoregressive parameter \eqn{\rho} less than 1, W is usually row-normalised
#' before implementing the SAR model. As in most cases, spatial weights matrix is very sparse, therefore W here should be
#' converted to a sparse matrix before imported into the `sar()` function to save computational burden and reduce computing
#' time. More specifically, W should be a column-oriented numeric sparse matrices of a `dgCMatrix` class defined in the
#' `Matrix` package. The converion between a dense numeric matrix and a sparse numeric matrix is made quite convenient through
#' the `Matrix` library.
#' @param Durbin `logical`. Estimate Durbin model (i.e. include spatial lags of `X` as predictors)? Default `FALSE`.
#' @param burnin The number of McMC samples to discard as the burnin period.
#' @param Nsim The total number of McMC samples to generate.
#' @param thinning MCMC thinning factor.
#' @param parameters.start A list with names "rho", "sigma2e", and "beta" corresponding to initial values for the model parameters
#' \eqn{\rho, \sigma^2_e} and the regression coefficients respectively.
#'
#' @return A `list`.
#' \describe{
#'  \item{cbetas}{A matrix with the MCMC samples of the draws for the coefficients.}
#'  \item{Mbetas}{A vector of estimated mean values of regression coefficients.}
#'  \item{SDbetas}{The standard deviations of estimated regression coefficients.}
#'  \item{Mrho}{The estimated mean of the lower-level spatial autoregressive parameter \eqn{\rho}.}
#'  \item{SDrho}{The standard deviation of the estimated lower-level spatial autoregressive parameter.}
#'  \item{Msigma2e}{The estimated mean of the lower-level variance parameter \eqn{\sigma^{2}_{e} }.}
#'  \item{SDsigma2e}{The standard deviation of the estimated lower-level variance parameter \eqn{\sigma^{2}_{e} }.}
#'  \item{DIC}{The deviance information criterion (DIC) of the fitted model.}
#'  \item{pd}{The effective number of parameters of the fitted model. }
#'  \item{Log_Likelihood}{The log-likelihood of the fitted model.}
#'  \item{R_Squared}{A pseudo R square model fit indicator.  }
#'  \item{impact_direct}{Summaries of the direct impact of a covariate effect on the outcome variable.  }
#'  \item{impact_idirect}{Summaries of the indirect impact of a covariate effect on the outcome variable. }
#'  \item{impact_total}{ Summaries of the total impact of a covariate effect on the outcome variable. }
#' }
#' @export
#'
#' @examples
#' data(landprice)
#' head(landprice)
#' data(land)
#'
#' # extract the land parcel level spatial weights matrix
#' library(spdep)
#' library(Matrix)
#' nb.25 <- spdep::dnearneigh(land,0,2500)
#' # to a weights matrix
#' dist.25 <- spdep::nbdists(nb.25,land)
#' dist.25 <- lapply(dist.25,function(x) exp(-0.5 * (x / 2500)^2))
#' mat.25 <- spdep::nb2mat(nb.25,glist=dist.25,style="W")
#' W <- as(mat.25,"dgCMatrix")
#'
#' ## run the sar() function
#' res.formula <- lnprice ~ lnarea + lndcbd + dsubway + dpark + dele +
#'                 popden + crimerate + as.factor(year)
#' betas= coef(lm(formula=res.formula,data=landprice))
#' pars=list(rho = 0.5, sigma2e = 2.0, betas = betas)
#' \donttest{
#' res <- sar(res.formula,data=landprice,W=W,
#'            burnin=500, Nsim=1000, thinning=1,
#'            parameters.start=pars)
#' summary(res)
#' }
sar <- function(formula, data = NULL, W, Durbin = FALSE,
                burnin = 5000, Nsim = 10000, thinning = 1,
                parameters.start = NULL) {

    ## check input data and formula
    frame <- check_formula(formula, data)
    X <- get_X_from_frame(frame)
    y <- get_y_from_frame(frame)

    if (any(is.na(y))) stop("NAs in dependent variable", call. = FALSE)
    if (any(is.na(X))) stop("NAs in independent variable", call. = FALSE)

    if(length(Durbin) != 1)           {stop("`Durbin` must be either TRUE or FALSE")}
    if(!(Durbin %in% c(TRUE, FALSE))) {stop("`Durbin` must be either TRUE or FALSE")}

    n <- nrow(X)

    check_matrix_dimensions(W,n,'Wrong dimensions for matrix W' )

    Xlabel <- colnames(X)
    Xnames <- setdiff(Xlabel, "(Intercept)")

    # 2025-05-20: attach lag X if `Durbin == TRUE`
    # TODO 2025-05-20: is it OK to treat factor dummies the same as any other column?
    lag_X <- if(Durbin) {W %*% X[, Xnames]}
    lag_X <- if(Durbin) {as.matrix(lag_X)}
    lag_X <- if(Durbin) {`colnames<-`(lag_X, paste0("lag_", Xnames))}

    X <- if(Durbin) {cbind(X, lag_X)} else {X}

    detval <- lndet_imrw(W)

    #start parameters
    if (! is.null(parameters.start)){
      if(is_there_parameter(parameters.start, "rho")) rho <- parameters.start$rho else rho<-0.5
      if(is_there_parameter(parameters.start, "sigma2e")) sigma2e<- parameters.start$sigma2e else sigma2e <-2.0
      if(is_there_parameter(parameters.start, "betas")) {
        betas <- parameters.start$betas
        if (dim(X)[2]!= length(betas) ) stop("Starting values for Betas have got wrong dimension", call. = FALSE)
      }
      else betas <- stats::coef(stats::lm(formula,data))
    }
    else{
      rho<-0.5
      sigma2e <-2.0
      betas <- if(Durbin) {stats::coef(stats::lm.fit(X, y))} else {stats::coef(stats::lm(formula,data))}
    }

    result <- sar_cpp_arma(X, y, W, detval, burnin, Nsim, thinning, rho, sigma2e, betas, Durbin )
      #.Call("HSAR_sar_cpp_arma", PACKAGE = 'HSAR', X, y, W, detval,
       #            burnin, Nsim, thinning, rho, sigma2e, betas )

    class(result) <- "mcmc_sar"

    result$cbetas<-put_labels_to_coefficients(result$cbetas, colnames(X))
    result$Mbetas<-put_labels_to_coefficients(result$Mbetas, colnames(X))
    result$SDbetas<-put_labels_to_coefficients(result$SDbetas, colnames(X))

    result$labels <- Xlabel
    result$call <-match.call()
    result$formula <- formula
    result$Durbin <- Durbin

    return(result)
}
