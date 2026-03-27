#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include "diagnostics.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Estimate the y hat
mat y_hat_sar(const mat& X, const mat& betas, double rho, const sp_mat& W ){

  mat Xb = X*betas.t();
  int n = X.n_rows;
  sp_mat I_sp = speye<sp_mat>(n,n);

  //mat invA = spsolve(A,I_nn, "lapack");
  sp_mat invA = I_sp+rho*W+pow(rho,2)*(W*W)+pow(rho,3)*(W*W*W);

  mat yhat = invA*Xb;

  return yhat;
}

// Log likelihood function of SAR
double SAR_loglikelihood(const mat& X, const mat& y, double rho,
                         const mat& betas, double sigma2e,
                         const mat& detval, const sp_mat& W){

  int n = X.n_rows;

  //find log( det(A) )
  uvec indt = find(detval.col(0) > rho, 1, "first");
  double logdetA = detval(indt[0],1);

  sp_mat I_sp = speye<sp_mat>(n,n);
  sp_mat A = I_sp - rho * W;

  mat Ay = A*y;

  mat Xb = X*trans(betas);
  mat AymXbAymXb= trans(Ay-Xb)*(Ay-Xb);

  double log_lik_sar = (-n/2)*(log(2*datum::pi)+log( pow(sigma2e,2)))
  + logdetA - AymXbAymXb(0,0)/(2*pow(sigma2e,2));

  return log_lik_sar;
}

// draw rho-values using griddy Gibbs and inversion
double SAR_draw_rho(const mat& detval, const mat& e0e0,
                    const mat& eded,const mat& e0ed,double sig){

    double rho_draw = 0;
    int nrho = detval.n_rows;

    vec rho_exist = detval.col(0);
    vec log_detrho = detval.col(1);

    vec iota = ones<vec>(nrho);

    // Calculate Log-S(e(rho*)) given rho*
    vec S_rho = e0e0(0,0)*iota + pow(rho_exist, 2) *eded(0,0) - 2*rho_exist*e0ed(0,0);

    // Calculate the Log-density
    vec log_den = log_detrho - (1.0/(2.0*sig)) * S_rho;
    double adj = log_den.max();
    vec t = rep(adj, nrho);
    log_den = log_den - t ;

    // the density
    vec den = exp(log_den);
    // the interval separating rho is h=0.001
    double h = 0.001;


    // Integration to calculate the normalized constant
    // using the  trapezoid rule
    double ISUM = h* ( sum(den)-den[0]/2.0-den[nrho-1]/2.0  );
    vec norm_den = (1.0/ISUM)*den;

    // cumulative density
    vec cumu_den = cumsum(norm_den);

    // Inverse sampling
    double rnd = runif(1,0,1)(0)*sum(norm_den);
    uvec ind = find(cumu_den <= rnd) ;
    int idraw = max(ind);

    if(idraw >= 0 && idraw < nrho)
        rho_draw = rho_exist(idraw);

    return rho_draw;
  }

// [[Rcpp::export]]

List sar_cpp_arma(arma::mat X, arma::vec y, arma::sp_mat W, arma::mat detval,
                  int burnin, int Nsim, int thinning,
                  float rho_start, float sigma2e_start, arma::vec betas_start,
                  bool Durbin){

  // if(Durbin) {throw std::runtime_error("here");}

  //Starting the MCMC SET UP

  int n = X.n_rows;
  int p = X.n_cols;

  //Prior distribution specifications
  //For Betas
  vec M0 = betas_start;

  mat T0 = 100.0 * eye<mat>(p,p);

  //completely non-informative priors
  double c0 = 0.01;
  double d0 = 0.01;
  // int c0(0.01);
  // int d0(0.01);

  //Store MCMC results
  mat Betas = zeros(Nsim-burnin, p);

  vec sigma2e = zeros(Nsim-burnin);

  vec rho = zeros(Nsim-burnin);

  //initial values for model parameters
  float sigma2e_i(sigma2e_start);
  float sigma2e_ip1(sigma2e_start);
  float rho_i(rho_start);
  float rho_ip1(rho_start);

  int ce( n/2 + c0 );

  //Fixed matrix manipulations during the MCMC loops

  mat XTX = trans(X) * X;
  mat invT0 = inv(T0);

  mat T0M0 = invT0 * M0;
  mat tX = trans(X);

  //some fixed values when updating rho
  mat beta0 = solve(X,y);

  mat e0 = y-X*beta0;
  mat e0e0 = trans(e0) * e0;

  vec Wy = W *y;
  mat betad = solve(X,Wy);
  mat ed = Wy - X * betad;

  mat eded = trans(ed) *ed;
  mat e0ed = trans(e0) *ed;

  // initialise A
  sp_mat I_sp = speye<sp_mat>(n,n);
  sp_mat A = I_sp - rho_i * W;

  // MCMC updating
  for(int i=1;i<=Nsim*thinning;i++) {
    // Gibbs sampler for updating Betas
    mat VV = (1.0/sigma2e_i) *XTX + invT0;
    // We use Cholesky decomption to inverse the covariance matrix
    mat vBetas = inv_sympd(VV) ;//chol2inv(chol(VV))

    // Define A=I-rho*W
    A = I_sp - rho_i * W;

    vec Ay = A*y;
    mat mBetas = vBetas * (tX * (1.0/sigma2e_i)*Ay+T0M0);

    // When the number of independent variables is large, drawing from multivariate
    // distribution could be time-comsuming
    mat cholV = trans(chol(vBetas));
    // draw form p-dimensitional independent norm distribution
    mat betas = rnorm(p);
    betas = mBetas + cholV * betas;

    //### Gibbs sampler for updating sigma2e
    mat Xb = X*betas;
    mat e = Ay - Xb;
    mat de = 0.5 * trans(e) * e + d0;

    sigma2e_ip1 = 1/Rf_rgamma(ce,1/de(0,0));
    rho_ip1 = SAR_draw_rho(detval, e0e0, eded, e0ed, sigma2e_ip1);

    rho_i = rho_ip1;
    sigma2e_i = sigma2e_ip1;

    if(i>(burnin*thinning))
    {
      if (i%thinning==0)
      {
        Betas.row(int(i-burnin*thinning)/thinning-1) = trans(betas);
        sigma2e[int(i-burnin*thinning)/thinning-1] = sigma2e_i;
        rho[int(i-burnin*thinning)/thinning-1] = rho_i;
      }
    }
  }

  // evaluate diagnostics
  vec log_lik_samples = zeros(Nsim-burnin);

  for(int i=0;i<(Nsim-burnin);i++)
  {
  log_lik_samples[i] = SAR_loglikelihood( X, y, rho[i], Betas.row(i),
                                    sigma2e[i], detval, W );
  }
  double log_lik_mean_theta = SAR_loglikelihood( X, y, mean( rho ),
                              mean( Betas ),
                             mean( sigma2e ), detval,
                             W );

  double dic, pd;
  diagnostic_dic_pd(log_lik_samples,log_lik_mean_theta, dic, pd);

  double r2 = diagnostic_Rsquared(y, y_hat_sar(X, mean( Betas ), mean( rho ),W ));

  mat direct, indirect, total;
  if(Durbin) {
    int half_p = (p - 1) / 2;

    mat means = mean( Betas );

    mat mean_betas = means.cols(1, half_p);
    mat mean_thetas = means.cols(half_p + 1, p - 1);

    diagnostic_impacts_Durbin(mean_betas, mean_thetas, mean(rho), W,
                              direct, indirect, total);
  } else {
    mat mean_betas = mean( Betas );

    diagnostic_impacts( mean_betas.cols(1, p - 1), mean( rho ),W , // 2025-05-27: don't include impact of `(Intercept)`
                        direct, indirect, total);
  }

  return List::create(Named("cbetas")= Betas,
                      Named("Mbetas")= mean( Betas ),
                      Named("SDbetas") = stddev( Betas ),
                      Named("crho") = rho,
                      Named("Mrho")= mean( rho ),
                      Named("SDrho") = stddev( rho ),
                      Named("csigma2e") = sigma2e,
                      Named("Msigma2e")= mean( sigma2e ),
                      Named("SDsigma2e") = stddev( sigma2e ),
                      Named("DIC") = dic,
                      Named("pD") = pd,
                      Named("Log_Likelihood") = log_lik_mean_theta,
                      Named("R_Squared") = r2,
                      Named("impact_direct") = direct,
                      Named("impact_indirect") = indirect,
                      Named("impact_total") = total,
                      Named("W") = W); // 2025-05-27: return W so we can take samples of impacts later

}



