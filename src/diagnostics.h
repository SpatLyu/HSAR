#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma;

void diagnostic_impacts(const mat& betas, double rho, const sp_mat& W,
                        mat& direct, mat& indirect, mat& total);

void diagnostic_impacts_Durbin(const mat& betas, const mat& thetas, double rho, const sp_mat& W,
                               mat& direct, mat& indirect, mat& total);

void diagnostic_dic_pd(vec log_likelihood_post_samples,
                       double log_likelihood_mean_theta,
                       double& DIC, double& pd);

double diagnostic_Rsquared(const mat& y, const mat& y_hat);

