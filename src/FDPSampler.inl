void FDPSampler::iteration_sample(std::mt19937 &rng){


	calculate_b();
	sample_cluster_labels(rng);
	update_weights(rng);

	draw_z(rng);
	X_fit << Z,X_K;
	V = (X_fit.transpose() * w.asDiagonal() * X_fit + correction_mat +  tau_matrix ).inverse();

	beta = V.llt().matrixL().toDenseMatrix() * sigma * z + V * X_fit.transpose() * w.asDiagonal() * y; 

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
	std::chi_squared_distribution<double> rchisq_tau(1 + (Q-P) );
	var = rchisq(rng);
	var = (s * (n - Q)) / var;
	sigma = sqrt(var);
	double temp_scale;
	temp_scale = (1 + beta.tail(Q-P).dot(beta.tail(Q-P)) / var ) / ((Q-P) + 1);
	temp_scale = sqrt(rchisq_tau(rng)) / ((1.0 + (Q-P)) * temp_scale) ;


	for(int k = 0; k< (Q-P); k += P_two){
		for(int p_ix =0; p_ix < P_two ; p_ix ++ ){
			tau_matrix(P + k+p_ix,P + k+p_ix) = temp_scale;
		}
	}

}

void FDPSampler::store_samples(Eigen::ArrayXXd &beta_samples,
						       Eigen::ArrayXd &sigma_samples,
							   Eigen::ArrayXXd &tau_samples,
							   Eigen::ArrayXXd &pi_samples,
							   Eigen::ArrayXd &alpha_samples,
							   Eigen::ArrayXXi &cluster_assignment){

	beta_samples.row(sample_ix) = beta.transpose().array();
	sigma_samples(sample_ix) = sigma;
	pi_samples.row(sample_ix) = pi.transpose();
	alpha_samples(sample_ix) = alpha;
	tau_samples(sample_ix,0) = pow(1 / tau_matrix(P,P),2) ;
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
	for(int k = 0; k < K; k++){
		b.col(k)  = log(pi(k))  -.5 * log(2 * M_PI * var) - 
			(.5  / var * pow(  (residual -  X  * beta.segment(start,P_two) ).array(),2)).array() - .5 * (beta.segment(start,P_two) * tau_matrix.block(start,start,P_two,P_two) * tau_matrix.block(start,start,P_two,P_two) *	(beta.segment(start,P_two)) ).array() ;
		start += P_two;
	}

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
		cluster_count(k) = (iter_cluster_assignment == k).count();
		if(cluster_count(k_)==0)
			correction_mat.diagonal().segment(P+k,P_two) = Eigen::VectorXd::Ones(P_two);
		X_K.block(0,k,n,P_two) = cluster_matrix.col(k_).asDiagonal() * X;
		k_ ++;
	}

}

void FDPSampler::update_weights(std::mt19937 &rng){

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
	for(int k = 0; k< (Q-P); k += P_two){
		for(int p_ix =0; p_ix < P_two ; p_ix ++ ){
			beta(P + k+p_ix) = rnorm(rng) + tau_0(p_ix);
			tau_matrix(P + k+p_ix,P + k+p_ix) = 1 / tau_0(p_ix);
		}
	}

}


