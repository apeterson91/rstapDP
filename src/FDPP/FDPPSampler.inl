void FDPPSampler::iteration_sample(std::mt19937 &rng){

	// update cluster labels
	calculate_b();
	log_message("b calculated");
	sample_cluster_labels(rng);
	log_message("labels sampled");
	update_weights(rng);
	log_message("weights updated");

	// create appropriately sized containers
	beta_temp.setZero(temp_Q);
	X_fit.setZero(n,temp_Q);
	X_fit.block(0,0,n,P) = Z;
	draw_z(rng);
	X_fit.block(0,P,n,temp_Q-P) = X_K;
	log_message("matrix designed");


	// calculate solution
	V = (X_fit.transpose() * w.asDiagonal() * X_fit + 
			nonzero_ics.transpose() * PenaltyMat * nonzero_ics ).inverse();

	log_message("matrix inverted");
	beta_temp = V.llt().matrixL().toDenseMatrix() * sigma * z + 
		V * X_fit.transpose() * w.asDiagonal() * y ;

	beta = nonzero_ics * beta_temp;

	adjust_beta(rng);
	log_message("beta sampled");


	// check for errors
	if(std::isnan(beta(0)) & flag ){
		Rcpp::Rcout << "There was an error inverting the design matrix" << std::endl;
		Rcpp::Rcout << "You may need to adjust the threshold parameter" << std::endl;
		Rcpp::Rcout << "Find out more in ?fdp_staplm.fit" << std::endl;
		flag = false;
	}

	// calculate predicted value
	yhat = X_fit * beta;
	log_message("yhat calculated");

	// draw variance terms
	draw_var(rng);
	log_message("sigma sampled");

}

void FDPPSampler::calculate_b(){

	b.setZero(n,K);
	residual = y - Z * beta.head(P);
	int start = P;
	for(int k = 0; k < K; k++){
		b.col(k)  =  log(pi(k)) - .5 * log(2 * M_PI ) - log(sigma * (1/w.array()))  -  
			(.5  *  precision * w.array() * pow( (residual -  X  * beta.segment(start,P_two) ).array(),2));
		start += P_two;
	}

}

void FDPPSampler::stick_break(std::mt19937 &rng){

	if(initializing){
		sftrabbit::beta_distribution<> rbeta(1,alpha);
		for(int k = 0; k < (K-1); k++){
			u(k) = rbeta(rng);
		}
		u(K-1) = 1;
	}else{
		for(int k = 0; k < (K-1); k++){
			u_posterior_beta_alpha(k) = 1 + cluster_count(k);
			u_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
			sftrabbit::beta_distribution<> rbeta(u_posterior_beta_alpha(k),u_posterior_beta_beta(k));
			u(k) = rbeta(rng);
		}
		u(K-1) = 1;
	}
	for(int k = 0; k < K; k++){
		pi(k) = k == 0 ?  u(k) : u(k) * (Eigen::ArrayXd::Ones(k) - u.head(k)).prod();
	}
}

void FDPPSampler::sample_cluster_labels(std::mt19937 &rng){

	cluster_matrix.setZero(n,K);

	for(int i = 0; i < n ; i ++){
		probs = b.row(i);
		probs = probs - probs.maxCoeff();
		for(int k = 0 ; k< K; k ++)
			probs(k) = probs(k) > log_factor ? exp(probs(k)) : 0;
		std::discrete_distribution<int> d(probs.data(),probs.data() + probs.size());
		iter_cluster_assignment(i) = d(rng);
		cluster_matrix(i,iter_cluster_assignment(i)) = 1.0;
	}

	log_message("cluster labels: assigned");
	adjust_zero_clusters(rng);


}
void FDPPSampler::adjust_zero_clusters(std::mt19937 &rng){

	std::gamma_distribution<double> rgamma(1,1);
	num_nonzero = K;

	for(int k = 0; k < K; k++){
		cluster_count(k) = (iter_cluster_assignment == k).count();
		if(cluster_count(k) <= threshold)
			num_nonzero --;
	}

	if(num_nonzero<=0)
		throw std::range_error("cannot have <= 0 nonzero clusters - error present");

	temp_Q = P + P_two * num_nonzero;
	nonzero_ics.setZero(Q,temp_Q);
	nonzero_ics.block(0,0,P,P) = Eigen::MatrixXd::Identity(P,P);
	X_K.setZero(n,P_two*num_nonzero);
	log_message("Finished preparing X_K Zero Matrix");

	int k_ = 0;
	for(int k = 0; k < K; k++){
		if(cluster_count(k)>threshold){
			X_K.block(0,k_*P_two,n,P_two) = cluster_matrix.col(k).asDiagonal() * X;
			nonzero_ics.block(P+k*P_two,P+k_*P_two,P_two,P_two) = Eigen::MatrixXd::Identity(P_two,P_two);
			k_ ++;
		}
	}
}

void FDPPSampler::update_weights(std::mt19937 &rng){

	stick_break(rng);
	if(!fix_alpha){
		posterior_b_alpha = 1.0 / alpha_b - (log(1 - u.array())).head(K-1).sum();
		posterior_b_alpha = 1.0 / posterior_b_alpha;
		std::gamma_distribution<double>  rgamma(posterior_a_alpha,posterior_b_alpha);
		alpha = rgamma(rng);
	}
}


void FDPPSampler::draw_var(std::mt19937 &rng){

	// scale for sigma precision
	residual = y - X_fit * beta_temp ;
	s = (residual.transpose() * w.asDiagonal()).dot(residual) * .5;
	s += .5 * (beta_temp.transpose() * nonzero_ics.transpose() * PenaltyMat * nonzero_ics).dot(beta_temp) ;
	std::gamma_distribution<double> rgamma(sigma_a + n * .5 + P_two*num_nonzero *.5, 1/( (1 / sigma_b) + s) );
	precision = rgamma(rng);
	sigma = sqrt(1 / precision);
	double temp_scale;
	PenaltyMat.setZero(Q,Q);
	for(int k = 0; k< K; k++){
		if(cluster_count(k)>threshold){
			temp_scale = calculate_penalty_scale(k,subset_one,false);
			std::gamma_distribution<double> rgamma_tau(tau_a + subset_one * .5 , 1/( (1/tau_b) + temp_scale) );
			unique_taus(k) = rgamma_tau(rng);
			update_penaltymat(k,subset_one,false);
			temp_scale = calculate_penalty_scale(k,subset_two,true);
			std::gamma_distribution<double> rgamma_tau2(tau_a + subset_two *.5 , 1/( (1/tau_b) + temp_scale) );
			unique_tau2s(k) = rgamma_tau2(rng);
			update_penaltymat(k,subset_two,true);
		}else{
			std::gamma_distribution<double> rgamma_tau_prior(tau_a,1/tau_b);
			unique_taus(k) = rgamma_tau_prior(rng);
			update_penaltymat(k,subset_one,false);
			std::gamma_distribution<double> rgamma_tau2_prior(tau_a,1/tau_b);
			unique_tau2s(k) = rgamma_tau_prior(rng);
			update_penaltymat(k,subset_two,true);
		}
	}
}

double FDPPSampler::calculate_penalty_scale(const int &k,const int &subset,const bool second){

	double out = 0;
	int diag_ix = P+k*P_two;
	diag_ix += second ? P_two - subset : 0;
	out = .5 * beta.segment(diag_ix,subset).squaredNorm();

	return(out);

}

void FDPPSampler::store_samples(Eigen::ArrayXXd &beta_samples,
						       Eigen::ArrayXd &sigma_samples,
							   Eigen::ArrayXXd &pi_samples,
							   Eigen::ArrayXXd &tau_samples,
							   Eigen::ArrayXXd &tau2_samples,
							   Eigen::ArrayXd &alpha_samples,
							   Eigen::ArrayXXi &cluster_assignment,
							   Eigen::ArrayXXd &yhat_samples){

	beta_samples.row(sample_ix) = beta.transpose().array();
	sigma_samples(sample_ix) = sigma;
	pi_samples.row(sample_ix) = pi.transpose();
	alpha_samples(sample_ix) = alpha;
	cluster_assignment.row(sample_ix) = iter_cluster_assignment;
	tau_samples.row(sample_ix) = unique_taus.transpose();
	tau2_samples.row(sample_ix) = unique_tau2s.transpose();
	yhat_samples.row(sample_ix) = yhat;
	sample_ix ++;
	for(int i = 0; i< n; i ++){
		for(int k = 0; k < K; k ++){
			if(cluster_matrix(i,k) == 1)
				P_matrix.row(i) = P_matrix.row(i) + cluster_matrix.col(k).transpose().array();
		}
	}

}

void FDPPSampler::initialize_beta(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	std::gamma_distribution<double> rgamma(tau_a,(1/tau_b));
	for(int p=0;p<Q;p++)
		beta(p) = rnorm(rng);

	for(int k = 0; k< K; k++){
		unique_taus(k) = rgamma(rng);
		update_penaltymat(k,subset_one,false);
		unique_tau2s(k) = rgamma(rng);
		update_penaltymat(k,subset_two,true);
	}
}

void FDPPSampler::draw_z(std::mt19937 &rng){

	z.setZero(temp_Q);
	std::normal_distribution<double> rnorm(0,1);
	for(int q = 0; q < temp_Q; q ++)
		z(q) = rnorm(rng);

}

void FDPPSampler::update_penaltymat(const int &k,const int &subset, const bool second){


	int diag_ix = P+k*P_two;
	diag_ix += second ? P_two - subset : 0;
	double penalty = second ? unique_tau2s(k) : unique_taus(k);
	PenaltyMat.block(diag_ix,diag_ix,subset,subset) = Eigen::MatrixXd::Identity(subset,subset) *  penalty;
}

void FDPPSampler::adjust_beta(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);

	for(int k = 0; k < K ; k++){
		if(cluster_count(k)<=threshold){
			for(int p = 0; p < P_two; p++)
				beta(P+P_two*k+p) = rnorm(rng)*sigma * sqrt(1/unique_taus(k));
		}
	}
}
