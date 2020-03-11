// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>
#include "FDPSampler.hpp"

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
Rcpp::List stapDP_backfit(const Eigen::VectorXd &y,
						  const Eigen::MatrixXd &Z,
						  const Eigen::VectorXd &X,
						  const double &tau_0,
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
	Eigen::MatrixXd P;
	P.setZero(y.rows(),y.rows());
	cluster_assignment.setZero(num_posterior_samples,y.rows());
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);

	const int chain = 1;
	int sample_ix = 0;

	FDPSampler sampler(y,Z,X,alpha_a,
					   alpha_b,tau_0,K,rng);


	sampler.iteration_sample(rng);

	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
		print_progress(iter_ix,burn_in,iter_max,chain);
		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
							      alpha_samples,cluster_assignment);
		}
	}

	for(int iter_ix = 0 ; iter_ix < num_posterior_samples; iter_ix ++){
		for(int i = 0; i < y.rows(); i ++){
			for(int j = 0; j < i; j ++)
				P(i,j) += cluster_assignment(iter_ix,i) == cluster_assignment(iter_ix,j) ? 1 : 0 ;
		}
	}

	P = P / num_posterior_samples;


    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = P);
}
