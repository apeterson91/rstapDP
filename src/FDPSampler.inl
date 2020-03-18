void FDPSampler::iteration_sample(std::mt19937 &rng){


	calculate_b();
	sample_cluster_labels(rng);
	update_weights(rng);

	draw_z(rng);
	X_fit << Z,X_K;
	V = (X_fit.transpose() * X_fit + tau_matrix ).inverse();

	beta = V * sigma * z + V * X_fit.transpose() * y; 
	if(isnan(beta(0))){
		Rcpp::Rcout << "pi: " << pi << std::endl;
		Rcpp::Rcout << "tau: " << (1/tau_matrix.diagonal().array()) << std::endl;
		Rcpp::Rcout << "sigma: " << sigma << std::endl;
		Rcpp::Rcout << "cluster assignment: " << cluster_matrix << std::endl;
	}


	draw_var(rng);
}

void FDPSampler::draw_z(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	for(int q = 0; q < Q; q ++)
		z(q) = rnorm(rng);

}

void FDPSampler::draw_var(std::mt19937 &rng){

	residual = y - X_fit * beta;
	s = residual.dot(residual) / (n - Q);
	std::chi_squared_distribution<double> rchisq(n - Q);
	var = rchisq(rng);
	var = (s * (n - Q)) / var;
	sigma = sqrt(var);
	int k_ = 0;
	double temp_scale;
	std::chi_squared_distribution<double> rchisq_tau(nu_0 + 1);

	for(int k = 0; k< (Q-P); k += P_two){
		k_ ++ ;
		for(int p_ix =0; p_ix < P_two ; p_ix ++ ){
			temp_scale = (nu_0  + pow(beta(P + k + p_ix),2) ) / (nu_0 + 1) ;
			tau_matrix(P + k + p_ix,P + k + p_ix) = 1 ;// 1.0 / (sqrt( (nu_0+1 * temp_scale  )  / rchisq_tau(rng) )); // invert the scaled inverse chi-square variate
		}
	}


}

void FDPSampler::store_samples(Eigen::ArrayXXd &beta_samples,
							  Eigen::ArrayXXd &tau_samples,
						       Eigen::ArrayXd &sigma_samples,
							   Eigen::ArrayXXd &pi_samples,
							   Eigen::ArrayXd &alpha_samples,
							   Eigen::ArrayXXi &cluster_assignment){

	beta_samples.row(sample_ix) = beta.transpose().array();
	tau_samples.row(sample_ix) = (1 / tau_matrix.diagonal().tail(Q-P).transpose().array());
	sigma_samples(sample_ix) = sigma;
	pi_samples.row(sample_ix) = pi.transpose();
	alpha_samples(sample_ix) = alpha;
	cluster_assignment.row(sample_ix) = iter_cluster_assignment;
	sample_ix ++;
	for(int i = 0; i< n; i ++){
		for(int j = 0; j < i ; j++)
			P_matrix(i,j) += iter_cluster_assignment(i) == iter_cluster_assignment(j) ? 1 : 0;
	}
}

void FDPSampler::stick_break(std::mt19937 &rng){

	if(initializing){
		sftrabbit::beta_distribution<> rbeta(1,alpha);
		for(int k = 0; k < (K-1); k++){
			u(k) = rbeta(rng);
		}
		u(K-1) = 1;
	}else{
		for(int k =0; k < (K-1); k++){
			sftrabbit::beta_distribution<> rbeta(u_posterior_beta_alpha(k),u_posterior_beta_beta(k));
			u(k) = rbeta(rng);
		}
		u(K-1) = 1;
	}
	for(int k = 0; k < K; k++)
		pi(k) = k == 0 ?  u(k) : u(k) * (Eigen::ArrayXd::Ones(k) - u.head(k)).prod();

}

void FDPSampler::calculate_b(){

	b.setZero(n,K);
	residual = y - Z * beta.head(P);
	int start = P;
	/*
	Rcpp::Rcout << " sigma: " << sigma << std::endl;
	Rcpp::Rcout << " pi: " << pi(0) << std::endl;
	Rcpp::Rcout << " delta: " << beta.segment(0,P) << std::endl;
	Rcpp::Rcout << " beta: " << beta.segment(P,P_two) << std::endl;
	Rcpp::Rcout << " tau_matrix \n" << tau_matrix << std::endl;
	Rcpp::Rcout << " tau_matrix transformed " << tau_matrix.diagonal().segment(start,P_two) << std::endl;
	*/
	for(int k = 0; k < K; k++){
		b.col(k)  = log(pi(k))   - (.5  / var * pow(  (residual -  X  * beta.segment(start,P_two) ).array(),2)).array() - 
			0.5 * ( pow(beta.segment(start,P_two).array(),2) * pow(tau_matrix.diagonal().segment(start,P_two).array(),2)   ).sum();//  -  
			//.5 * P_two *  log(2*M_PI) -  P_two * (log( (1.0 / tau_matrix.diagonal().segment(start,P_two).array() ) )).sum() -  
	//		(nu_0 *  (.5 * tau_matrix.diagonal().segment(start,P_two).array() * 
	//				  tau_matrix.diagonal().segment(start,P_two).array() )  - (1 + nu_0 /2 ) * log(tau_matrix.diagonal().segment(start,P_two).array() ) ).sum();  
		start += P_two;
	}
	//Rcpp::Rcout << " b(0,0):" <<  b(0,0) << std::endl;

}

void FDPSampler::sample_cluster_labels(std::mt19937 &rng){

	cluster_matrix.setZero(n,K);
	double log_factor = log(pow(10,-16)) - log(n);
	for(int i = 0; i < n ; i ++){
		probs = b.row(i);
		probs = probs - probs.maxCoeff();
		for(int k = 0 ; k< K; k ++)
			probs(k) = probs(k) > log_factor ? exp(probs(k)) : 0;
		std::discrete_distribution<int> d(probs.data(),probs.data() + probs.size());
		iter_cluster_assignment(i) = d(rng);
		cluster_matrix(i,iter_cluster_assignment(i)) = 1.0;
	}
	int k_ = 0;
	for(int k = 0; k < (Q-P); k+= P_two){
		X_K.block(0,k,n,P_two) = cluster_matrix.col(k_).asDiagonal() * X;
		k_ ++;
	}

}

void FDPSampler::update_weights(std::mt19937 &rng){

	for(int k = 0; k < K; k++)
		cluster_count(k) = (iter_cluster_assignment == k).count();

	for(int k = 0; k < K; k++){
		u_posterior_beta_alpha(k) = 1 + cluster_count(k);
		u_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
	}

	stick_break(rng);
	posterior_b_alpha = 1.0 / alpha_b - (log(1 - u.array())).head(K-1).sum();
	std::gamma_distribution<double>  rgamma(posterior_a_alpha,posterior_b_alpha);
	alpha = rgamma(rng);

}

void FDPSampler::initialize_beta(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	std::chi_squared_distribution<double> rchisq(nu_0);
	for(int k = 0; k< (Q-P); k += P_two){
		for(int p_ix =0; p_ix < P_two ; p_ix ++ ){
			tau_matrix(P + k+p_ix,P + k+p_ix) =  sqrt( rchisq(rng) / nu_0) ; // invert the inverted chi square variate
			beta(P + k+p_ix) = rnorm(rng)* (1.0 / tau_matrix(P + k + p_ix, P + k + p_ix)); 
		}
	}

}


