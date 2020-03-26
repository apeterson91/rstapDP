// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>
#include "FDPSampler.hpp"
#include "FDPPSampler.hpp"

void print_progress(const int &iter_ix, const int &warm_up, const int &iter_max, const int &chain){
  if(iter_max > 20){
      if((iter_ix) % (int)round(.1 * iter_max) == 0 || iter_ix == 1 || iter_ix == (warm_up + 1) ){
          int progress = (int)round(iter_ix * 100 / iter_max);
          std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
          Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
      }
  }
  else{
          int progress = (int)round(iter_ix * 100 / iter_max);
          std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
          Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
  }

}

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
//
// [[Rcpp::export]]
Rcpp::List stapDP_fit(const Eigen::VectorXd &y,
						  const Eigen::MatrixXd &Z,
						  const Eigen::MatrixXd &X,
						  const Eigen::ArrayXd tau_0,
						  const double &alpha_a,
						  const double &alpha_b,
						  const int &K,
						  const int &iter_max,
						  const int &burn_in,
						  const int &thin,
						  const int &seed,
						  const int &num_posterior_samples
						  ) {

    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);
	
	// create sample containers
	Eigen::ArrayXXd beta_samples;
	Eigen::ArrayXd sigma_samples;
	Eigen::ArrayXd alpha_samples;
	Eigen::ArrayXXd pi_samples;
	Eigen::ArrayXXi cluster_assignment;
	cluster_assignment.setZero(num_posterior_samples,y.rows());
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + X.cols()*K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);

	const int chain = 1;
	int sample_ix = 0;

	FDPSampler sampler(y,Z,X,alpha_a,
					   alpha_b,tau_0,K,
					   rng);


	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
		print_progress(iter_ix,burn_in,iter_max,chain);
		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
							      alpha_samples,cluster_assignment);
		}
	}



    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = sampler.P_matrix / num_posterior_samples );
}

//' Penalized Functional Dirichlet Process Regression
//'
//' @param y vector of outcomes 
//' @param Z matrix of subject covariates 
//' @param X matrix of distance spline basis 
// [[Rcpp::export]]
Rcpp::List stappDP_fit(const Eigen::VectorXd &y,
						  const Eigen::MatrixXd &Z,
						  const Eigen::MatrixXd &X,
						  const double &nu_0,
						  const double &alpha_a,
						  const double &alpha_b,
						  const int &K,
						  const int &iter_max,
						  const int &burn_in,
						  const int &thin,
						  const int &seed,
						  const int &num_posterior_samples
						  ) {

    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);
	
	// create sample containers
	Eigen::ArrayXXd beta_samples;
	Eigen::ArrayXd sigma_samples;
	Eigen::ArrayXd alpha_samples;
	Eigen::ArrayXXd pi_samples;
	Eigen::ArrayXXd tau_samples;
	Eigen::ArrayXXi cluster_assignment;
	cluster_assignment.setZero(num_posterior_samples,y.rows());
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + X.cols()*K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);
	tau_samples.setZero(num_posterior_samples,X.cols()*K);

	const int chain = 1;
	int sample_ix = 0;

	FDPPSampler sampler(y,Z,X,alpha_a,
					   alpha_b,nu_0,K,rng);


	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
		print_progress(iter_ix,burn_in,iter_max,chain);
		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
								  tau_samples,alpha_samples,
								  cluster_assignment);
		}
	}



    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("tau") = tau_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = sampler.P_matrix / num_posterior_samples );
}
