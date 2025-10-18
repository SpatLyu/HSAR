#include <iostream>
#include <RcppArmadillo.h>
#include "diagnostics.h"
#include "gibbs_method.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// Estimate the y hat
mat y_hat_hsar(const mat& X, const mat& betas,
               double rho, const sp_mat& W,
               const sp_mat& Z,const mat& Mus){

  mat Xb = X*betas.t();
  int n = X.n_rows;

  sp_mat I_sp = speye<sp_mat>(n,n);

  //mat invA = spsolve(A,I_nn, "lapack");
  sp_mat invA = I_sp+rho*W+pow(rho,2)*(W*W)+pow(rho,3)*(W*W*W);

  mat Zg = Z* Mus.t();

  mat yhat = invA*(Xb+Zg);

  return yhat;
}

// Log likelihood function of HSAR
double HSAR_loglikelihood(const mat& X, const mat& y, double rho,
                          const mat& betas, const mat& us,
                          const vec& Unum, int Utotal, double sigma2e,
                          const mat& detval, const sp_mat& W){

  int n = X.n_rows;

  //find log( det(A) )
  uvec indt = find(detval.col(0) > rho, 1, "first");
  double logdetA = detval(indt[0],1);

  sp_mat I_sp = speye<sp_mat>(n,n);
  sp_mat A = I_sp - rho * W;

  mat Ay = A*y;

  mat Xb = X*trans(betas);

  mat Zu;
  mat nus = trans(us);
  for(int j=0;j<Utotal;j++){
    mat temp_add = repmat(nus.row(j),Unum[j], 1 );

    Zu.insert_rows( Zu.n_rows, temp_add);
  }

  mat crosprod_AymXbmZu= trans(Ay-Xb-Zu)*(Ay-Xb-Zu);

  double log_lik_sar = (-n/2)*(log(2*datum::pi)+log( pow(sigma2e,2)))
  + logdetA - crosprod_AymXbmZu(0,0)/(2*pow(sigma2e,2));

  return log_lik_sar;
}

// [[Rcpp::export]]

List hsar_cpp_arma(arma::mat X, arma::vec y, arma::sp_mat W, arma::sp_mat M,
                   arma::sp_mat Z, arma::mat detval, arma::mat detvalM,
                   arma::vec Unum, int burnin, int Nsim, int thinning,
                   float rho_start, float lambda_start,
                   float sigma2e_start, float sigma2u_start, arma::vec betas_start,
                   bool Durbin){

  //Starting the MCMC SET UP
  //arma_rng::set_seed(124);

  int n = X.n_rows;
  int p = X.n_cols;
  int Utotal = M.n_rows;

  //Prior distribution specifications
  //For Betas
  vec M0 = betas_start;

  mat T0 = 100.0 * eye<mat>(p,p);

  //completely non-informative priors
  double c0 = 0.01;
  double d0 = 0.01;
  double a0 = 0.01;
  double b0 = 0.01;
  // int c0(0.01);
  // int d0(0.01);
  // int a0(0.01);
  // int b0(0.01);

  //Store MCMC results
  mat Betas = zeros(Nsim-burnin, p);
  mat Us = zeros(Nsim-burnin, Utotal);

  vec sigma2e = zeros(Nsim-burnin);
  vec sigma2u = zeros(Nsim-burnin);

  vec rho = zeros(Nsim-burnin);
  vec lambda = zeros(Nsim-burnin);

  //initial values for model parameters
  float sigma2e_i(sigma2e_start), sigma2e_ip1(sigma2e_start);
  float sigma2u_i(sigma2u_start), sigma2u_ip1(sigma2u_start);

  float rho_i(rho_start), rho_ip1(rho_start);
  float lambda_i(lambda_start), lambda_ip1(lambda_start);

  int ce(n/2 + c0);
  int au(Utotal/2 + a0);

  //Fixed matrix manipulations during the MCMC loops
  mat XTX = trans(X) * X;
  mat invT0 = inv(T0);

  mat T0M0 = invT0 * M0;
  mat tX = trans(X);
  mat Zfull(Z);
  mat tZ = trans(Zfull);

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
  sp_mat I_sp_A = speye<sp_mat>(n,n);
  sp_mat A = I_sp_A - rho_i * W;

 //  initialise B
  sp_mat I_sp_B = speye<sp_mat>(Utotal,Utotal);
  sp_mat B = I_sp_B - lambda_i * M;

  // MCMC updating
  mat us = zeros(Utotal);
  for(int i=1;i<=Nsim*thinning;i++) {
    // Gibbs sampler for updating Betas
    mat VV = (1.0/sigma2e_i) *XTX + invT0;
    // We use Cholesky decomption to inverse the covariance matrix
    mat vBetas = inv_sympd(VV) ;//chol2inv(chol(VV))

    //clock_t A_start_time = clock();
    // Define A=I-rho*W
    A = I_sp_A - rho_i * W;

    vec Ay = A*y;
    vec ZUs = Z*us ;

    mat mBetas = vBetas * (tX * (1.0/sigma2e_i)*(Ay-ZUs)+T0M0);

    // When the number of independent variables is large, drawing from multivariate
    // distribution could be time-comsuming
    mat cholV = trans(chol(vBetas));
    // draw form p-dimensitional independent norm distribution
    mat betas = Rcpp::rnorm(p);
    betas = mBetas + cholV * betas;

    // Gibbs sampler for updating U. Now, us are spatially dependent.

    // Define and update B=I-lambda*M
    B = I_sp_B - lambda_i * M;
    mat vU = tZ * (1.0/sigma2e_i)*Z+ trans(B) * (1.0/sigma2u_i)*B;
    vU = inv_sympd(vU); //vU <- chol2inv(chol(vU))

    vec Xb = X*betas;
    mat mU = vU * (tZ * (1.0/sigma2e_i)*(Ay-Xb ));

    // When the number of higher level units is large, drawing from multivariate
    // distribution could be time-comsuming
    cholV = trans(chol(vU));

    // draw form J-dimensitional independent norm distribution
    us = Rcpp::rnorm(Utotal);
    us = mU + cholV * us;

    // Gibbs sampler for updating sigma2e
    mat Zu;
    for(int j=0;j<Utotal;j++){
      mat temp_add = repmat(us.row(j),Unum[j], 1 );

      Zu.insert_rows( Zu.n_rows, temp_add);
    }

    mat e = Ay - Zu -Xb;
    mat de = 0.5 * trans(e) * e + d0;

    sigma2e_ip1 = 1/Rf_rgamma(ce,1/de(0,0));

   // Gibbs sampler for updating sigma2u
   vec Bus = B*us;
   mat bu = 0.5 * trans(Bus) * Bus + b0;
   sigma2u_ip1 = 1/Rf_rgamma(au,1/bu(0,0));

   // Giddy Gibbs integration and inverse sampling for rho
    mat betau = solve(X,Zu);
    mat eu = Zu-X*betau;

    mat eueu = trans(eu) * eu;
    mat e0eu = trans(e0) * eu;
    mat edeu = trans(ed) * eu;

    //clock_t rho_start_time = clock();
    rho_ip1 = HSAR_draw_rho(detval, e0e0, eded, eueu, e0ed, e0eu, edeu, sigma2e_ip1);

    // Giddy Gibbs integration and inverse sampling for lambda

    mat uu = trans(us) *us;
    mat uMu = trans(us) * M * us;
    mat Mu = M * us;
    mat uMMu = trans(Mu) * Mu;

    lambda_ip1 = HSAR_draw_lambda(detvalM, uu, uMu, uMMu, sigma2u_ip1);

    rho_i = rho_ip1;lambda_i = lambda_ip1;
    sigma2e_i = sigma2e_ip1;sigma2u_i = sigma2u_ip1;

    if(i>(burnin*thinning))
    {
      if (i%thinning==0)
      {
        Betas.row(int(i-burnin*thinning)/thinning-1) = trans(betas);
        Us.row(int(i-burnin*thinning)/thinning-1) = trans(us);

        sigma2e[int(i-burnin*thinning)/thinning-1] = sigma2e_i;
        sigma2u[int(i-burnin*thinning)/thinning-1] = sigma2u_i;

        rho[int(i-burnin*thinning)/thinning-1] = rho_i;
        lambda[int(i-burnin*thinning)/thinning-1] = lambda_i;
      }
    }
  }

  // Diagnostics
  vec log_lik_samples = zeros(Nsim-burnin);

  for(int i=0;i<(Nsim-burnin);i++)
  {
  log_lik_samples[i] = HSAR_loglikelihood(X, y, rho[i], Betas.row(i),
                                          Us.row(i),Unum,Utotal,
                                          sigma2e[i], detval, W );
  }
  double log_lik_mean_theta = HSAR_loglikelihood(X, y, mean( rho ),
                                                 mean(Betas),mean(Us), Unum,Utotal,
                                                 mean(sigma2e), detval, W);

  double dic, pd;
  diagnostic_dic_pd(log_lik_samples,log_lik_mean_theta, dic, pd);

  double r2 = diagnostic_Rsquared(y, y_hat_hsar(X, mean( Betas ), mean( rho ), W, Z ,mean( Us ) ));

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
                      Named("crho") = rho, // 2025-05-06: why not return the draws of all the things
                      Named("Mrho")= mean( rho ),
                      Named("SDrho") = stddev( rho ),
                      Named("clambda") = lambda,
                      Named("Mlambda")= mean( lambda ),
                      Named("SDlambda") = stddev( lambda ),
                      Named("csigma2e")= sigma2e,
                      Named("Msigma2e")= mean( sigma2e ),
                      Named("SDsigma2e") = stddev( sigma2e ),
                      Named("csigma2u")= mean( sigma2u ),
                      Named("Msigma2u")= mean( sigma2u ),
                      Named("SDsigma2u") = stddev( sigma2u ),
                      Named("cus")=  Us,
                      Named("Mus")= mean( Us ),
                      Named("SDus") = stddev( Us ),
                      Named("DIC") = dic,
                      Named("pD") = pd,
                      Named("Log_Likelihood") = log_lik_mean_theta,
                      Named("R_Squared") = r2,
                      Named("impact_direct") = direct,
                      Named("impact_indirect") = indirect,
                      Named("impact_total") = total,
                      Named("W") = W); // 2025-05-27: return W so we can take samples of impacts later

}


