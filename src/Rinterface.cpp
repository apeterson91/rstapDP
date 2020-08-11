// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>
#include "FDPP/FDPPSampler.hpp"


#include "auxiliary/print_function.hpp"

// [[Rcpp::depends(RcppEigen)]]

//' Penalized Functional Dirichlet Process Linear Regression
//'
//' @param y a vector of continuous outcomes
//' @param Z a matrix of population level confounders
//' @param X a matrix of spatial temporal aggregated predictors
//' @param S penalty matrix for stap parameters
//' @param w a vector of weights for weighted regression
//' @param alpha_a alpha gamma prior shape hyperparameter
//' @param alpha_b alpha gamma prior scale hyperparameter
//' @param sigma_a precision gamma prior shape hyperparameter
//' @param sigma_b precision gamma prior scale hyperparameter
//' @param tau_a penalty gamma prior shape hyperparameter
//' @param tau_b penalty gamma prior scale hyperparameter
//' @param K truncation number
//' @param num_penalties number of penalty matrices accounted for in S
//' @param iter_max maximum number of iterations
//' @param burn_in number of burn in iterations
//' @param thin number by which to thin samples
//' @param seed rng initializer
//' @param num_posterior_samples total number of posterior samples
//' @param fix_alpha  boolean value that determines whether or not to fix alpha in sampler
// [[Rcpp::export]]
Rcpp::List stappDP_fit(const Eigen::VectorXd &y,
					   const Eigen::MatrixXd &Z,
					   const Eigen::MatrixXd &X,
					   const Eigen::MatrixXd &S,
					   const Eigen::VectorXd &w,
					   const double &alpha_a,
					   const double &alpha_b,
					   const double &sigma_a,
					   const double &sigma_b,
					   const double &tau_a,
					   const double &tau_b,
					   const int &K,
					   const int &num_penalties,
					   const int &iter_max,
					   const int &burn_in,
					   const int &thin,
					   const int &seed,
					   const int &num_posterior_samples,
					   const bool &fix_alpha
					   ) {

    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);
	
	// create sample containers
	Eigen::ArrayXXd beta_samples;
	Eigen::ArrayXd sigma_samples;
	Eigen::ArrayXd alpha_samples;
	Eigen::ArrayXXd pi_samples;
	Eigen::ArrayXXd yhat_samples;
	Eigen::ArrayXXd tau_samples;
	Eigen::ArrayXXi cluster_assignment;
	cluster_assignment.setZero(num_posterior_samples,y.rows());
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + X.cols()*K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);
	tau_samples.setZero(num_posterior_samples,K*num_penalties);
	yhat_samples.setZero(num_posterior_samples,y.rows());

	const int chain = 1;

	FDPPSampler sampler(y,Z,X,S,w, alpha_a,alpha_b,tau_a,tau_b,
						sigma_a,sigma_b,K,num_penalties,fix_alpha,rng);


	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){

		print_progress(iter_ix,burn_in,iter_max,chain);

		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
								  tau_samples,alpha_samples,
								  cluster_assignment,yhat_samples);
		}
	}



    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("tau") = tau_samples,
							  Rcpp::Named("yhat") = yhat_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = sampler.P_matrix / num_posterior_samples );
}

typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> SpMat;
#include "FDPP_mer/FDPPSampler_mer.hpp"

//' Penalized Functional Dirichlet Process Linear Mixed Effects Regression
//'
//' @param y a vector of continuous outcomes
//' @param Z a matrix of population level confounders
//' @param X a matrix of spatial temporal aggregated predictors
//' @param W a design matrix for group specific terms
//' @param S penalty matrix for stap parameters
//' @param w a vector of weights for weighted regression
//' @param alpha_a alpha gamma prior shape hyperparameter
//' @param alpha_b alpha gamma prior scale hyperparameter
//' @param sigma_a precision gamma prior shape hyperparameter
//' @param sigma_b precision gamma prior scale hyperparameter
//' @param tau_a penalty gamma prior shape hyperparameter
//' @param tau_b penalty gamma prior scale hyperparameter
//' @param K truncation number
//' @param num_penalties number of penalty matrices accounted for in S
//' @param iter_max maximum number of iterations
//' @param burn_in number of burn in iterations
//' @param thin number by which to thin samples
//' @param seed rng initializer
//' @param num_posterior_samples total number of posterior samples
//' @param fix_alpha  boolean value that determines whether or not to fix alpha in sampler
// [[Rcpp::export]]
Rcpp::List stappDP_mer_fit(const Eigen::VectorXd &y,
						   const Eigen::MatrixXd &Z,
						   const Eigen::MatrixXd &X,
						   const Eigen::ArrayXXd &W,
						   const Eigen::MatrixXd &S,
						   const Eigen::VectorXd &w,
						   const SEXP &subj_mat_,
						   const Eigen::ArrayXi &subj_n,
						   const double &alpha_a,
						   const double &alpha_b,
						   const double &sigma_a,
						   const double &sigma_b,
						   const double &tau_a,
						   const double &tau_b,
						   const int &K,
						   const int &num_penalties,
						   const int &iter_max,
						   const int &burn_in,
						   const int &thin,
						   const int &seed,
						   const int &num_posterior_samples,
						   const bool &fix_alpha
							)
{

    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);

	//Map SparseMatrix
	const SpMat subj_mat(Rcpp::as<MSpMat>(subj_mat_));
	
	// create sample containers
	Eigen::ArrayXXd beta_samples;
	Eigen::ArrayXd sigma_samples;
	Eigen::ArrayXd alpha_samples;
	Eigen::ArrayXXd pi_samples;
	Eigen::ArrayXXd yhat_samples;
	Eigen::ArrayXXd tau_samples;
	Eigen::ArrayXXi cluster_assignment;
	Eigen::ArrayXXd b_samples;
	Eigen::ArrayXXd D_samples;
	cluster_assignment.setZero(num_posterior_samples,subj_mat.cols());
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + X.cols()*K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);
	tau_samples.setZero(num_posterior_samples,K*num_penalties);
	yhat_samples.setZero(num_posterior_samples,y.rows());
	b_samples.setZero(num_posterior_samples,W.cols()*subj_mat.cols());
	D_samples.setZero(num_posterior_samples,W.cols()*W.cols());


	const int chain = 1;

	FDPPSampler_mer sampler(y,Z,X,W,S,w,
							subj_mat,subj_n,
							alpha_a,
							alpha_b,tau_a,tau_b,
							sigma_a,sigma_b,K,
							num_penalties,fix_alpha,rng);

	sampler.iteration_sample(rng);

	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){

		print_progress(iter_ix,burn_in,iter_max,chain);

		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
								  tau_samples,alpha_samples,
								  cluster_assignment,yhat_samples,
								  b_samples,D_samples);
		}
	}



    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("tau") = tau_samples,
							  Rcpp::Named("yhat") = yhat_samples,
							  Rcpp::Named("subj_b") = b_samples,
							  Rcpp::Named("subj_D") = D_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = sampler.P_matrix / num_posterior_samples );
}
