#' @title Hierarchical SAR model estimation
#' @description
#'  The specification of a HSAR model is as follows:
#'  \deqn{y_{i,j} = \rho *\mathbf{W}_i *\mathbf{y} + \mathbf{x}^\prime_{i,j} * \mathbf{\beta} +
#'    \mathbf{z}^\prime_j * \mathbf{\gamma} + \theta_j + \epsilon_{i,j}  }{y_i,j = \rho * W_i * y + x'_i,j * \beta +z'_j * \gamma + \theta_j + \epsilon_i,j}
#'  \deqn{\theta_j = \lambda * \mathbf{M}_j * \mathbf{\theta} + \mu_j  }{\theta_j = \lambda * M_j * \theta + \mu_j}
#'  \deqn{\epsilon_{i,j} \sim N(0,\sigma_e^2), \hspace{2cm} \mu_j \sim N(0,\sigma_u^2)}{\epsilon_i,j ~ N(0,\sigma^2_e); \mu_j ~ N(0,\sigma^2_u)}
#'  where \eqn{i=1,2,...,n_j} and \eqn{j=1,2,...,J} are indicators of lower- and higher-level spatial units. \eqn{n_j} is the number of lower-level units in the \eqn{j-th} higher level unit and \eqn{\sum_{j=1}^J=\mathbf{N}}{Sum(j=1,J)=N}. \eqn{\mathbf{x}^\prime_{i,j}}{x'_i,j} and \eqn{\mathbf{z}^\prime_j}{z'_j} represent vectors of lower- and higher-level independent variables. \eqn{\mathbf{\beta}}{\beta} and \eqn{\mathbf{\gamma}}{\gamma} are regression coefficients to estimate. \eqn{\mathbf{\theta}}{\theta}, a \eqn{N \times J}{N * J} vector of higher-level random effects, also follows a simultaneous autoregressive process. \eqn{\mathbf{W}}{W} and \eqn{\mathbf{M}}{M} are two spatial weights matrices (or neighbourhood connection matrices) at the lower and higher levels, defining how spatial units at each level are connected. \eqn{\rho} and \eqn{\lambda} are two spatial autoregressive parameters measuring the strength of the dependencies/correlations at the two spatial scales.

#'  A succinct matrix formulation of the model is,
#'  \deqn{\mathbf{y} = \rho * \mathbf{W} * \mathbf{y} + \mathbf{X} * \mathbf{\beta} +
#'    \mathbf{Z} * \mathbf{\gamma} + \Delta * \mathbf{\theta} + \mathbf{\epsilon} }{y = \rho * W * y + X * \beta + Z * \gamma + \Delta * \theta + \epsilon}

#'  \deqn{\mathbf{\theta} = \lambda * \mathbf{M} * \mathbf{\theta} + \mathbf{\mu}}{\theta = \lambda * M * \theta + \mu }

#'  It is also useful to note that the HSAR model nests a standard (random intercept) multilevel model model when \eqn{\rho} and \eqn{\lambda} are both equal to zero and a standard spaital econometric model when \eqn{\lambda} and \eqn{\sigma^2_u} are both equal to zero.
#'
#' @note
#' In order to use the hsar() function, users need to specify the two spatial weights matrices W and M and the random effect design matrix \eqn{\delta}. However, it is very easy to extract such spatial weights matrices from spatial data using the package \pkg{spdep}. Geographic distance-based or contiguity-based spatial weights matrix for both spatial points data and spatial polygons data are available in the \pkg{spdep} package.
#' Before the extraction of W and M, it is better to first sort the data using the higher-level unit identifier. Then, the random effect design matrix can be extracted simply (see the following example) and so are the two spatial weights matrices. Make sure the order of higher-level units in the weights matrix M is in line with that in the \eqn{\delta} matrix.
#' Two simpler versions of the HSAR model can also be fitted using the hsar() function. The first is a HSAR model with \eqn{\lambda} equal to zero, indicating an assumption of independence in the higher-level random effect \eqn{\mathbf{\theta}}. The second is a HSAR with \eqn{\rho} equal to zero, indicating an independence assumption in the outcome variable conditioning on the hgiher-level random effect. This model is useful in situations where we are interested in the neighbourhood/contextual effect on individual's outcomes and have good reasons to suspect the effect from geographical contexts upon individuals to be dependent. Meanwhile we have no information on how lower-level units are connnected.
#'
#' @references
#' Dong, G. and Harris, R. 2015. Spatial Autoregressive Models for Geographically Hierarchical Data Structures. \emph{Geographical Analysis}, 47:173-191.
#'
#' LeSage, J. P., and R. K. Pace. (2009). \emph{Introduction to Spatial Econometrics}. Boca Raton, FL: CRC Press/Taylor & Francis.
#'
#' @param formula A symbolic description of the model to fit. A formula for the covariate part of the model using the syntax of the lm() function fitting standard linear regression models. Neither the response variable nor the explanatory variables are allowed to contain NA values.
#' @param data A `data.frame` containing variables used in the formula object.
#' @param W The N by N lower-level spatial weights matrix or neighbourhood matrix where N is the total number of lower-level spatial units. The formulation of W could be based on geographical distances separating units or based on geographical contiguity. To ensure the maximum value of the spatial autoregressive parameter \eqn{\rho} less than 1,
#' W should be row-normalised before running the HSAR model. As in most cases, spatial weights matrix is very sparse, therefore W here should be converted to a sparse matrix before imported into the `hsar()` function to save computational burden and reduce computing time. More specifically, W should be a column-oriented numeric sparse matrices of
#' a `dgCMatrix` class defined in the `Matrix` package. The converion between a dense numeric matrix and a sparse numeric matrix is made quite convenient through the `Matrix`library.
#' @param M The J by J higher-level spatial weights matrix or neighbourhood matrix where J is the total number of higher-level spatial units. Similar with W, the formulation of M could be based on geographical distances separating units or based on geographical contiguity. To ensure the maximum value of the spatial autoregressive parameter \eqn{\lambda}
#' less than 1, M is also row-normalised before running the HSAR model. As with W, M should also be a column-oriented numeric sparse matrices.
#' @param Delta The N by J random effect design matrix that links the J by 1 higher-level random effect vector back to the N by 1 response variable under investigation. It is simply how lower-level units are grouped into each high-level units with columns of the matrix being each higher-level units. As with W and M, \eqn{\delta} should also be a
#' column-oriented numeric sparse matrices.
#' @param Durbin `logical`. Estimate Durbin model (i.e. include spatial lags of `X` as predictors)? Default `FALSE`.
#' @param burnin The number of MCMC samples to discard as the burnin period.
#' @param Nsim The total number of MCMC samples to generate.
#' @param thinning MCMC thinning factor.
#' @param parameters.start A list with names "rho", "lambda", "sigma2e", "sigma2u" and "beta" corresponding to initial values for the model parameters \eqn{\rho, \lambda, \sigma^2_e, \sigma^2_u} and the regression coefficients respectively.
#'
#' @return A `list`.
#' \describe{
#'  \item{cbetas}{A matrix with the MCMC samples of the draws for the coefficients.}
#'  \item{Mbetas}{A vector of estimated mean values of regression coefficients. }
#'  \item{SDbetas}{The standard deviations of estimated regression coefficients.}
#'  \item{crho}{A vector with the MCMC samples of the draws for the lower-level spatial autoregressive parameter.}
#'  \item{Mrho}{The estimated mean of the lower-level spatial autoregressive parameter \eqn{\rho}.}
#'  \item{SDrho}{The standard deviation of the estimated lower-level spatial autoregressive parameter.}
#'  \item{clambda}{A vector with the MCMC samples of the draws for the higher-level spatial autoregressive parameter.}
#'  \item{Mlambda}{The estimated mean of the higher-level spatial autoregressive parameter \eqn{\lambda}.}
#'  \item{SDlambda}{The standard deviation of the estimated higher-level spatial autoregressive parameter.}
#'  \item{csigma2e}{A vector with the MCMC samples of the draws for the lower-level variance parameter.}
#'  \item{Msigma2e}{The estimated mean of the lower-level variance parameter \eqn{\sigma^2_e}.}
#'  \item{SDsigma2e}{The standard deviation of the estimated lower-level variance parameter \eqn{\sigma^{2}_{e} }.}
#'  \item{csigma2u}{A vector with the MCMC samples of the draws for the higher-level variance parameter.}
#'  \item{Msigma2u}{The estimated mean of the higher-level variance parameter \eqn{\sigma^2_u}.}
#'  \item{SDsigma2u}{The standard deviation of the estimated higher-level variance parameter \eqn{\sigma^2_u}.}
#'  \item{cus}{A matrix with the MCMC samples of the draws of \eqn{\theta}.}
#'  \item{Mus}{Mean values of \eqn{\theta} }
#'  \item{SDus}{Standard deviation of \eqn{\theta} }
#'  \item{DIC}{The deviance information criterion (DIC) of the fitted model.}
#'  \item{pd}{The effective number of parameters of the fitted model.  }
#'  \item{Log_Likelihood}{The log-likelihood of the fitted model. }
#'  \item{R_Squared}{A pseudo R square model fit indicator.   }
#'  \item{impact_direct}{Summaries of the direct impact of a covariate effect on the outcome variable.}
#'  \item{impact_idirect}{Summaries of the indirect impact of a covariate effect on the outcome variable. }
#'  \item{impact_total}{Summaries of the total impact of a covariate effect on the outcome variable.}
#' }
#' @export
#'
#' @examples
#' library(spdep)
#'
#' # Running the hsar() function using the Beijing land price data
#' data(landprice)
#'
#' # load shapefiles of Beijing districts and land parcels
#' data(Beijingdistricts)
#' data(land)
#'
#' plot(Beijingdistricts,border="green")
#' plot(land,add=TRUE,col="red",pch=16,cex=0.8)
#'
#' # Define the random effect matrix
#' model.data <- landprice[order(landprice$district.id),]
#' head(model.data,50)
#'
#' # the number of individuals within each neighbourhood
#' MM <- as.data.frame(table(model.data$district.id))
#' # the total number of neighbourhood, 100
#' Utotal <- dim(MM)[1]
#' Unum <- MM[,2]
#' Uid <- rep(c(1:Utotal),Unum)
#'
#' n <- nrow(model.data)
#' Delta <- matrix(0,nrow=n,ncol=Utotal)
#' for(i in 1:Utotal) {
#'   Delta[Uid==i,i] <- 1
#' }
#' rm(i)
#' # Delta[1:50,1:10]
#' Delta <- as(Delta,"dgCMatrix")
#'
#' # extract the district level spatial weights matrix using the queen's rule
#' nb.list <- spdep::poly2nb(Beijingdistricts)
#' mat.list <- spdep::nb2mat(nb.list,style="W")
#' M <- as(mat.list,"dgCMatrix")
#'
#' # extract the land parcel level spatial weights matrix
#' nb.25 <- spdep::dnearneigh(land,0,2500)
#' # to a weights matrix
#' dist.25 <- spdep::nbdists(nb.25,land)
#' dist.25 <- lapply(dist.25,function(x) exp(-0.5 * (x / 2500)^2))
#' mat.25 <- spdep::nb2mat(nb.25,glist=dist.25,style="W")
#' W <- as(mat.25,"dgCMatrix")
#'
#'
#' ## run the hsar() function
#' res.formula <- lnprice ~ lnarea + lndcbd + dsubway + dpark + dele +
#'   popden + crimerate + as.factor(year)
#'
#' betas= coef(lm(formula=res.formula,data=landprice))
#' pars=list( rho = 0.5,lambda = 0.5, sigma2e = 2.0, sigma2u = 2.0, betas = betas )
#'
#' \donttest{
#'   res <- hsar(res.formula, data=landprice, W=W, M=M, Delta=Delta,
#'               burnin=500, Nsim=1000, thinning = 1, parameters.start=pars)
#'   summary(res)
#'
#'   # visualise the district level random effect
#'   groups <- sdsfun::discretize_vector(res$Mus,n = 4,method = "natural")
#'   palette <- RColorBrewer::brewer.pal(4, "Blues")
#'   plot(Beijingdistricts,col=palette[groups],border="grey")
#' }
hsar <- function(formula, data = NULL, W=NULL, M=NULL, Delta, Durbin = FALSE,
                 burnin=5000, Nsim=10000, thinning=1, parameters.start = NULL) {

    ## check input data and formula
    frame <- check_formula(formula, data)
    X <- get_X_from_frame(frame)
    y <- get_y_from_frame(frame)

    if (any(is.na(y))) stop("NAs in dependent variable", call. = FALSE)
    if (any(is.na(X))) stop("NAs in independent variable", call. = FALSE)

    if( is.null(W) & is.null(M) ) stop('Both weight matrices can not be NULL', call. = FALSE)

    n <- nrow(Delta)
    p <- ncol(Delta)

    if( !is.null(W) ) check_matrix_dimensions(W,n,'Wrong dimensions for matrix W' )
    if( !is.null(M) ) check_matrix_dimensions(M,p,'Wrong dimensions for matrix M' )

    if(length(Durbin) != 1)           {stop("`Durbin` must be either TRUE or FALSE")}
    if(!(Durbin %in% c(TRUE, FALSE))) {stop("`Durbin` must be either TRUE or FALSE")}

    Xlabel <- colnames(X)
    Xnames <- setdiff(Xlabel, "(Intercept)")

    # 2025-05-20: attach lag X if `Durbin == TRUE`
    # TODO 2025-05-20: is it OK to treat factor dummies the same as any other column?
    lag_X <- if(Durbin) {W %*% X[, Xnames]}
    lag_X <- if(Durbin) {as.matrix(lag_X)}
    lag_X <- if(Durbin) {`colnames<-`(lag_X, paste0("lag_", Xnames))}

    X <- if(Durbin) {cbind(X, lag_X)} else {X}

    Unum <- apply(Delta,2,sum) # TODO 2025-05-20: eventually allow for multiple membership somehow

    #start parameters
    if (! is.null(parameters.start)){
      if(is_there_parameter(parameters.start, "rho")) rho <- parameters.start$rho else rho<-0.5
      if(is_there_parameter(parameters.start, "sigma2e")) sigma2e<- parameters.start$sigma2e else sigma2e <-2.0
      if(is_there_parameter(parameters.start, "rho")) lambda <- parameters.start$lambda else lambda<-0.5
      if(is_there_parameter(parameters.start, "sigma2e")) sigma2u<- parameters.start$sigma2u else sigma2u <-2.0
      if(is_there_parameter(parameters.start, "betas")) {
        betas <- parameters.start$betas
        if (dim(X)[2]!= length(betas) ) stop("Starting values for Betas have got wrong dimension", call. = FALSE)
      }
      else betas <- stats::coef(stats::lm.fit(X, y)) # 2025-07-15: fixed this long ago in `sar.R` but not here -GA
    }
    else{
      rho <- 0.5
      lambda <- 0.5
      sigma2e <- 2.0
      sigma2u <- 2.0
      # betas <- stats::coef(stats::lm(formula,data)) # 2025-07-03: omits betas for lag X vars if Durbin = TRUE -GA
      betas <- stats::lm.fit(X, y)$coefficients # 2025-07-03: I thought I tried this before, though -GA # 2025-07-15: see above -GA
    }

    ## Call various models
    # Special case where rho =0 ; dependent regional effect
    if (is.null(W)){
      detval <- lndet_imrw(M)
      result <- hsar_cpp_arma_rho_0(X, y, M, Delta, detval, Unum, burnin, Nsim, thinning, lambda, sigma2e, sigma2u, betas)
        #.Call("HSAR_hsar_cpp_arma_rho_0", PACKAGE = 'HSAR', X, y, M, Delta, detval, Unum,
         #             burnin, Nsim, thinning, lambda, sigma2e, sigma2u, betas)
      class(result) <- "mcmc_hsar_rho_0"
    }
    # Special case where lamda =0 ; independent regional effect
    if ( is.null(M)){
      detval <- lndet_imrw(W)
      result <- hsar_cpp_arma_lambda_0(X, y, W, Delta, detval, Unum, burnin, Nsim, thinning, rho, sigma2e, sigma2u, betas, Durbin)
        #.Call("HSAR_hsar_cpp_arma_lambda_0", PACKAGE = 'HSAR', X, y, W, Delta, detval, Unum,
         #             burnin, Nsim, thinning, rho, sigma2e, sigma2u, betas)
      class(result) <- "mcmc_hsar_lambda_0"
    }
    # Full HSAR model
    if ( (!is.null(M)) & (!is.null(W))){
      detval <- lndet_imrw(W)
      detvalM <- lndet_imrw(M)
      result <- hsar_cpp_arma(X, y, W, M, Delta, detval, detvalM, Unum, burnin, Nsim, thinning, rho, lambda, sigma2e, sigma2u, betas, Durbin)
        #.Call("HSAR_hsar_cpp_arma", PACKAGE = 'HSAR', X, y, W, M, Delta, detval, detvalM, Unum,
         #             burnin, Nsim, thinning, rho, lambda, sigma2e, sigma2u, betas)
      class(result) <- "mcmc_hsar"
    }

    result$cbetas<-put_labels_to_coefficients(result$cbetas, colnames(X))
    result$Mbetas<-put_labels_to_coefficients(result$Mbetas, colnames(X))
    result$SDbetas<-put_labels_to_coefficients(result$SDbetas, colnames(X))

    if(!is.null(W)) { # 2025-05-27: name the impacts here and avoid complications arising from `Durbin == TRUE` downstream
      result$impact_direct   <- put_labels_to_coefficients(result$impact_direct, Xnames)
      result$impact_indirect <- put_labels_to_coefficients(result$impact_direct, Xnames)
      result$impact_total    <- put_labels_to_coefficients(result$impact_direct, Xnames)
    }

    result$labels <- Xlabel
    result$call <- match.call()
    result$formula <- formula
    result$Durbin <- Durbin

    return(result)
}
