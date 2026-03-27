#include "diagnostics.h"

// Calculate the summary measures of the impacts (using aproximation for inv(SW) )
// TODO 2025-05-27: pass in SW as an argument to make mcmc_impacts faster, I guess
// TODO 2025-05-27: consider increasing number of powers of rho * W for accuracy
void diagnostic_impacts(const mat& betas, double rho, const sp_mat& W,
                        mat& direct, mat& indirect, mat& total) {

  int n = W.n_rows;

  sp_mat I_sp = speye<sp_mat>(n,n);

  sp_mat SW = I_sp+rho*W+pow(rho,2)*(W*W)+pow(rho,3)*(W*W*W)+pow(rho,4)*(W*W*W*W)+pow(rho,5)*(W*W*W*W*W);

  vec d( SW.diag() );

  direct = sum(d)/n * betas;
  total = accu(SW)/n * betas;
  indirect = total - direct;
}

// 2025-05-22: as a non-C++ programmer I don't usually think about scope this way -GA
void diagnostic_impacts_Durbin(const mat& betas, const mat& thetas, double rho, const sp_mat& W,
                               mat& direct, mat& indirect, mat& total) {

  int n = W.n_rows;

  sp_mat I_sp = speye<sp_mat>(n,n);

  sp_mat SW = I_sp+rho*W+pow(rho,2)*(W*W)+pow(rho,3)*(W*W*W)+pow(rho,4)*(W*W*W*W)+pow(rho,5)*(W*W*W*W*W);

  vec d( SW.diag() );

  direct = sum(d)/n * betas;
  total = accu(SW)/n * betas;
  // indirect = total - direct;

  SW = W * SW;

  vec t( SW.diag() );

  direct += sum(t)/n * thetas;
  total += accu(SW)/n * thetas;
  indirect = total - direct;
}


// Estimate Deviance Information Criterion (DIC)
void diagnostic_dic_pd(vec log_likelihood_post_samples,double log_likelihood_mean_theta, double& DIC, double& pd) {

  double D_hat = -2*log_likelihood_mean_theta;

  double D_bar = -2*mean(log_likelihood_post_samples);

  pd = D_bar - D_hat;
  DIC = D_bar + pd ;
}

// Estimate R squared
double diagnostic_Rsquared(const mat& y, const mat& y_hat){

  double y_mean = mean( y.col(0) ) ;
  double yhat_mean = mean( y_hat.col(0) );

  int p = y.n_rows;

  double r2 = sum(pow(y_hat.col(0)-yhat_mean*ones<vec>(p),2)) / sum(pow(y.col(0)-y_mean*ones<vec>(p),2));

  return r2;
}
