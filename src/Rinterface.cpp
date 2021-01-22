// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>
#include "FDPP/FDPPSampler.hpp"
#include <progress.hpp>
 

#include "auxiliary/print_function.hpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppProgress)]]

//' Penalized Functional Dirichlet Process Linear Regression with N observations
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
//' @param chain chain label
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
					   const int chain,
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
//' fits a functional dirichlet process linear mixed effects regression model
//' with N observations and n subjects
//'
//' @param y a vector of continuous outcomes
//' @param Z a matrix of population level confounders
//' @param X a matrix of spatial temporal aggregated predictors
//' @param W a design matrix for group specific terms
//' @param w a vector of weights for weighted regression
//' @param subj_mat_ N x n sparse matrix used to aggregate subject observations
//' @param subj_n n x 1 vector of integers representing how many observations correspond to each subject
//' @param alpha_a alpha gamma prior shape hyperparameter
//' @param alpha_b alpha gamma prior scale hyperparameter
//' @param sigma_a precision gamma prior shape hyperparameter
//' @param sigma_b precision gamma prior scale hyperparameter
//' @param tau_a penalty gamma prior shape hyperparameter
//' @param tau_b penalty gamma prior scale hyperparameter
//' @param K truncation number
//' @param iter_max maximum number of iterations
//' @param burn_in number of burn in iterations
//' @param thin number by which to thin samples
//' @param seed rng initializer
//' @param chain chain label
//' @param num_posterior_samples total number of posterior samples
//' @param fix_alpha  boolean value that determines whether or not to fix alpha in sampler
//' @param logging boolean parameter indicating whether or not a single iteration should be run with print messages indicating successful completion of the Sampler's sub modules
//' @param summarize_yhat boolean value indicating whether a single mean vector of yhat values should be returned instead of a N X num samples matrix. Useful in situations where N is large.
// [[Rcpp::export]]
Rcpp::List stappDP_mer_fit(const Eigen::VectorXd &y,
						   const Eigen::MatrixXd &Z,
						   const Eigen::MatrixXd &X,
						   const Eigen::ArrayXXd &W,
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
						   const int &iter_max,
						   const int &burn_in,
						   const int &thin,
						   const int &seed,
						   const int &chain,
						   const int &num_posterior_samples,
						   const bool &fix_alpha,
						   const bool &logging,
						   const bool &summarize_yhat
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
	tau_samples.setZero(num_posterior_samples,K*X.cols());
	yhat_samples.setZero(num_posterior_samples,y.rows());
	b_samples.setZero(num_posterior_samples,W.cols()*subj_mat.cols());
	D_samples.setZero(num_posterior_samples,W.cols()*W.cols());

	Progress p(0,false);
	FDPPSampler_mer sampler(y,Z,X,W,w,
							subj_mat,subj_n,
							alpha_a,
							alpha_b,tau_a,tau_b,
							sigma_a,sigma_b,K,
							fix_alpha,logging,rng);
	if(logging){
		sampler.iteration_sample(rng);
		return(Rcpp::List::create(Rcpp::Named("log status") = 0));
	}

	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){

		print_progress(iter_ix,burn_in,iter_max,chain);

		sampler.iteration_sample(rng);
		if(iter_ix % 100 == 0){
			if(Progress::check_abort())
				return(Rcpp::List::create(Rcpp::Named("log status") = 1));
		}

		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
								  tau_samples,alpha_samples,
								  cluster_assignment,yhat_samples,
								  b_samples,D_samples);
		}
	}

	if(summarize_yhat){
		yhat_samples = yhat_samples.rowwise().mean();
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


#include "FDPP_mer/FDPPSamplerdecomp.hpp"
//' Penalized Functional Dirichlet Process Linear Mixed Effects Regression with Between-Within Decomposition
//'
//' fits a functional dirichlet process linear mixed effects regression model
//' with N observations and n subjects using a between-within effects decomposition on the STAP-DP design matrix 
//'
//' @param y a vector of continuous outcomes
//' @param Z a matrix of population level confounders
//' @param X_b Matrix of between subject spatial temporal aggregated predictor covariates
//' @param X_w Matrix of within subject spatial temporal aggregated predictor covariates
//' @param W a design matrix for group specific terms
//' @param S_b penalty matrix corresponding to between subject covariate matrix
//' @param S_w penalty matrix corresponding to within subject covariate matrix
//' @param w a vector of weights for weighted regression
//' @param subj_mat_ N x n sparse matrix used to aggregate subject observations
//' @param subj_n n x 1 vector of integers representing how many observations correspond to each subject
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
//' @param chain chain label
//' @param fix_alpha  boolean value that determines whether or not to fix alpha in sampler
// [[Rcpp::export]]
Rcpp::List stappDP_merdecomp(const Eigen::VectorXd &y,
						   const Eigen::MatrixXd &Z,
						   const Eigen::MatrixXd &X_b,
						   const Eigen::MatrixXd &X_w,
						   const Eigen::ArrayXXd &W,
						   const Eigen::MatrixXd &S_b,
						   const Eigen::MatrixXd &S_w,
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
						   const int &chain,
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
	Eigen::ArrayXXd tau_samples_b;
	Eigen::ArrayXXd tau_samples_w;
	Eigen::ArrayXXi cluster_assignment;
	Eigen::ArrayXXd b_samples;
	Eigen::ArrayXXd D_samples;
	cluster_assignment.setZero(num_posterior_samples,subj_mat.cols());
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + (X_w.cols() + X_b.cols())*K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);
	tau_samples_b.setZero(num_posterior_samples,K*num_penalties);
	tau_samples_w.setZero(num_posterior_samples,K*num_penalties);
	yhat_samples.setZero(num_posterior_samples,y.rows());
	b_samples.setZero(num_posterior_samples,W.cols()*subj_mat.cols());
	D_samples.setZero(num_posterior_samples,W.cols()*W.cols());

	FDPPSamplerdecomp sampler(y,Z,X_b,W,S_b,w,
							subj_mat,subj_n,
							alpha_a,
							alpha_b,tau_a,tau_b,
							sigma_a,sigma_b,K,
							num_penalties,
							fix_alpha,rng,
							X_w,S_w);

	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){

		print_progress(iter_ix,burn_in,iter_max,chain);

		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,pi_samples,
								  tau_samples_b,tau_samples_w,
								  alpha_samples,cluster_assignment,
								  yhat_samples,b_samples,D_samples);
		}
	}



    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("tau_b") = tau_samples_b,
							  Rcpp::Named("tau_w") = tau_samples_w,
							  Rcpp::Named("yhat") = yhat_samples,
							  Rcpp::Named("subj_b") = b_samples,
							  Rcpp::Named("subj_D") = D_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = sampler.P_matrix / num_posterior_samples );
}



#include "VI_DP/VI.hpp"

//'  Spatial Temporal Aggregated Predictor Functional Dirichlet Process Regression via Variational Inference
//' 
//' @export
//' @param y vector of regression outcomes
//' @param X matrix of regression covariates
//' @param tau_a penalty precision gamma shape hyperparameter  
//' @param tau_b penalty precision gamma scale hyperparameter  
//' @param sigma_a residual precision gamma shape hyperparameter  
//' @param sigma_b residual precision gamma scale hyperparameter  
//' @param max_iter maximum number of iterations
//' @param num_samples number of samples to draw from approximate posterior
//' @param seed random number generator seed initializer
// [[Rcpp::export]]
Rcpp::List VI_lm(const Eigen::VectorXd &y,
				 const Eigen::MatrixXd &X,
				 const double &tau_a,
				 const double &tau_b,
				 const double &sigma_a,
				 const double &sigma_b,
				 const int &max_iter,
				 const int &num_samples,
				 const int seed){

    // set seed
    std::mt19937 rng;
    rng = std::mt19937(seed);
	Eigen::ArrayXd bounds;
	Eigen::ArrayXd sigma_samples;
	Eigen::ArrayXd tau_samples;
	Eigen::ArrayXXd yhat_samples;
	Eigen::ArrayXXd beta_samples;
	//regression constants,predictions, coefficients
	const int Q = X.cols();
	const int N = X.rows();

	// containers
	bounds.setZero(max_iter);
	sigma_samples.setZero(num_samples);
	tau_samples.setZero(num_samples);
	yhat_samples.setZero(num_samples,N);
	beta_samples.setZero(num_samples,Q);

	VI sampler(y,X,tau_a,tau_b,sigma_a,sigma_b,max_iter,rng);
	  
	while(sampler.has_not_converged(bounds)){
		sampler.descend_gradient();
	}
	sampler.draw_pars(rng,num_samples,beta_samples,sigma_samples,tau_samples,yhat_samples);



	return Rcpp::List::create(Rcpp::Named("Bound") = bounds,
								Rcpp::Named("beta") = beta_samples,
								Rcpp::Named("sigma") = sigma_samples,
								Rcpp::Named("tau") = tau_samples);
}


// RcppParallel barebones example
/*
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

struct Test : public Worker {
  tbb::enumerable_thread_specific<bool> printonce;
  Test() : printonce(false) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    tbb::enumerable_thread_specific<bool>::reference p = printonce.local();
    if(!p) { // print once per thread
      std::cout << 1;
      p= true;
    }
  }
};

// [[Rcpp::export(rng = false)]]
void test() {
  Test x{};
  parallelFor(0, 10000, x);
}
*/
