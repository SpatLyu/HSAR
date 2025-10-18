#' Estimate impact per MCMC draw of coefficients
#' @param object an `mcmc_hsar` or `mcmc_hsar_lambda_0` or `mcmc_sar` object
#' @param p `numeric`. Maximum weight matrix power. Default `5`.
#' @return a list of three matrices named `direct`, `indirect`, and `total` giving impact estimates for each MCMC draw
#' @export
mcmc_impacts <- function(object, p = 5) {
  betas <- object$cbetas
  betas <- as.matrix(betas)
  betas <- betas[, setdiff(colnames(betas), "(Intercept)")]

  if(object$Durbin) {
    thetas <- betas[,  grepl("^lag_", colnames(betas)), drop = FALSE]
    betas  <- betas[, !grepl("^lag_", colnames(betas)), drop = FALSE]
    thetas <- `colnames<-`(thetas, stringr::str_remove(colnames(thetas), "^lag_"))
  }

  # take first p powers of wates
  wates <- object$W
  n     <- nrow(wates)
  eyes  <- Matrix::Diagonal(n)
  mates <- list(eyes)
  w8pow <- eyes
  i     <- 0
  while(i < p) {
    w8pow <- w8pow %*% wates
    mates <- c(mates, list(w8pow))
    i <- i + 1
  }

  # first p powers of each sample of rho, divided by n
  rhopo <- outer(as.numeric(object$crho), 0:p, `^`) / n

  diags <- lapply(mates, Matrix::diag)
  disum <- sapply(diags, sum)
  dirho <- rhopo %*% disum
  direc <- betas * as.numeric(dirho)

  # m8sum <- sapply(mates, sum) # these are all equal tho?
  m8sum <- rep(n, p + 1)
  m8rho <- rhopo %*% m8sum
  total <- betas * as.numeric(m8rho)

  # also impacts of lag vars
  if(object$Durbin) {
    mates <- lapply(mates, `%*%`, y = wates)

    liags <- lapply(mates, Matrix::diag)
    lisum <- sapply(liags, sum)
    lirho <- rhopo %*% lisum
    lirec <- thetas * as.numeric(lirho)
    direc <- direc + lirec

    lotal <- thetas * as.numeric(m8rho)
    total <- total + lotal
  }


  list(direct = direc, indirect = total - direc, total = total)
}
