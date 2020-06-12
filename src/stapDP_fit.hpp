
//
// via the exports attribute we tell Rcpp to make this function
// available from R

//' Functional Dirichlet Process Linear Regression
//' @param y a vector of continuous outcomes
//' @param Z a matrix of population level confounders
//' @param X a matrix of spatial temporal aggregated predictors
//' @param w  vector of weights for weighted regression
//' @param tau_0 prior variance for STAP parameters
//' @param alpha_a alpha gamma prior hyperparameter
//' @param alpha_b alpha gamma prior hyperparameter
//' @param K truncation number
//' @param iter_max maximum number of iterations
//' @param burn_in number of burn in iterations
//' @param thin number by which to thin samples
//' @param seed rng initializer
//' @param num_posterior_samples number of final samples
// [[Rcpp::export]]
Rcpp::List stapDP_fit(const Eigen::VectorXd &y,
						  const Eigen::MatrixXd &Z,
						  const Eigen::MatrixXd &X,
						  const Eigen::VectorXd &w,
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
	Eigen::ArrayXXd tau_samples;
	Eigen::ArrayXXd pi_samples;
	Eigen::ArrayXXi cluster_assignment;
	cluster_assignment.setZero(num_posterior_samples,y.rows());
	tau_samples.setZero(num_posterior_samples,1);
	alpha_samples.setZero(num_posterior_samples);
	beta_samples.setZero(num_posterior_samples,Z.cols() + X.cols()*K);
	sigma_samples.setZero(num_posterior_samples);
	pi_samples.setZero(num_posterior_samples,K);

	const int chain = 1;
	int sample_ix = 0;

	FDPSampler sampler(y,Z,X,w,alpha_a,
					   alpha_b,tau_0,K,
					   rng);


	for(int iter_ix = 1; iter_ix <= iter_max; iter_ix ++){
		print_progress(iter_ix,burn_in,iter_max,chain);
		sampler.iteration_sample(rng);
		if(iter_ix > burn_in && (iter_ix % thin == 0)){
			sampler.store_samples(beta_samples,sigma_samples,tau_samples,
								  pi_samples,alpha_samples,cluster_assignment);
		}
	}



    return Rcpp::List::create(Rcpp::Named("beta") = beta_samples,
							  Rcpp::Named("pi") = pi_samples,
							  Rcpp::Named("tau") = tau_samples,
                              Rcpp::Named("sigma") = sigma_samples,
							  Rcpp::Named("alpha") = alpha_samples,
							  Rcpp::Named("cluster_assignment") = cluster_assignment,
							  Rcpp::Named("PairwiseProbabilityMat") = sampler.P_matrix / num_posterior_samples );
}
