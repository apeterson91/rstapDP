void FDPPSampler::iteration_sample(std::mt19937 &rng){

	calculate_b();
	sample_cluster_labels(rng);
	update_weights(rng);

	draw_z(rng);
	X_fit << Z,X_K;
	V = (X_fit.transpose() * w.asDiagonal() * X_fit + tau_var_matrix + S ).inverse();

	beta = V.llt().matrixL().toDenseMatrix() * z * sigma + V * X_fit.transpose() * w.asDiagonal() *  y; 
	if(std::isnan(beta(0)) & flag ){
		Rcpp::Rcout << "things are NaN" << std::endl;
		Rcpp::Rcout << " V block" << V.block(0,0,5,5) << std::endl;
		Rcpp::Rcout << "w.diagonal" << w.head(5)  << std::endl;
		flag = false;
	}
	yhat = X_fit * beta;
	/*
	if( beta.squaredNorm()>10000){
		Rcpp::Rcout << " large beta: " << std::endl;
		Rcpp::Rcout << "cluster count" << cluster_count.transpose() << std::endl;
		Rcpp::Rcout << "beta" << beta.tail(Q-P).transpose() << std::endl;
		Rcpp::Rcout << "V" << V.block(P,P,P_two*K,P_two*K)   << std::endl;
		Rcpp::Rcout << "Lz sigma" << V.llt().matrixL().toDenseMatrix() * z * sigma << std::endl;
		Rcpp::Rcout << "VX^T" << V  * X_fit.transpose() * w.asDiagonal() * y << std::endl;
		Rcpp::Rcout << "taus: " <<  unique_taus_sq.transpose() << std::endl;
		Rcpp::Rcout << "z: " <<  z << std::endl;
		Rcpp::Rcout << "L: " <<  V.llt().matrixL().toDenseMatrix() << std::endl;
		Rcpp::Rcout << "sigma: " <<  sigma << std::endl;
	}
	*/
	/*
	if((cluster_count==0).count()>0){
		Rcpp::Rcout << "Zero member Cluster:" << std::endl;
		Rcpp::Rcout << cluster_count << std::endl;
		Rcpp::Rcout << "beta" << beta.tail(Q-P).transpose() << std::endl;
		Rcpp::Rcout << "unique_taus" << unique_taus_sq.sqrt() << std::endl;
	}
	*/

	draw_var(rng);

}

void FDPPSampler::stick_break(std::mt19937 &rng){

	if(initializing){
		sftrabbit::beta_distribution<> rbeta(1,alpha);
		for(int k = 0; k < (K-1); k++){
			u(k) = rbeta(rng);
		}
		u(K-1) = 1;
	}else{
		for(int k =0; k < (K-1); k++){
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

void FDPPSampler::calculate_b(){

	b.setZero(n,K);
	residual = y - Z * beta.head(P);
	int start = P;
	for(int k = 0; k < K; k++){
		b.col(k)  =  log(pi(k)) - .5 * log(2 * M_PI * var) -  
			(.5  / var * pow( (residual -  X  * beta.segment(start,P_two) ).array(),2));
		start += P_two;
	}

}

void FDPPSampler::sample_cluster_labels(std::mt19937 &rng){

	cluster_matrix.setZero(n,K);
	int cntr = 0;
	for(int i = 0; i < n ; i ++){
		probs = b.row(i);
		probs = probs - probs.maxCoeff();
		for(int k = 0 ; k< K; k ++)
			probs(k) = probs(k) > log_factor ? exp(probs(k)) : 0;
		std::discrete_distribution<int> d(probs.data(),probs.data() + probs.size());
		iter_cluster_assignment(i) = d(rng);
		cluster_matrix(i,iter_cluster_assignment(i)) = 1.0;
		cntr++;
	}
	int k_ = 0;
	std::chi_squared_distribution<double> rchisq(nu_0);
	for(int k = 0; k < (Q-P); k+= P_two){
		cluster_count(k_) = (iter_cluster_assignment == k_).count();
		X_K.block(0,k,n,P_two) = cluster_matrix.col(k_).asDiagonal() * X;
		if(cluster_count(k_)==0){
			unique_taus_sq(k_) = nu_0 / rchisq(rng);
			tau_var_matrix.block(P+P_two*k_,P+P_two*k_,P_two,P_two) = Eigen::MatrixXd::Identity(P_two,P_two) * (1 / unique_taus_sq(k_));
		}
		k_ ++;
	}

}

void FDPPSampler::update_weights(std::mt19937 &rng){

	stick_break(rng);
	posterior_b_alpha = 1.0 / alpha_b - (log(1 - u.array())).head(K-1).sum();
	posterior_b_alpha = 1.0 / posterior_b_alpha;
	std::gamma_distribution<double>  rgamma(posterior_a_alpha,posterior_b_alpha);
	alpha = rgamma(rng);

}

void FDPPSampler::draw_z(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	for(int q = 0; q < Q; q ++)
		z(q) = rnorm(rng);

}

void FDPPSampler::draw_var(std::mt19937 &rng){

	residual = y - X_fit * beta;
	s = ((residual.transpose() * w.asDiagonal() * residual) + beta.segment(P,Q-P).dot( tau_var_matrix.block(P,P,Q-P,Q-P) *  beta.segment(P,Q-P)) ) / (w.sum() + Q - P);
	std::chi_squared_distribution<double> rchisq(n + Q-P);
	std::chi_squared_distribution<double> rchisq_tau(nu_0 + P_two);
	var = rchisq(rng);
	var = (s * (n+Q-P)) / var;
	sigma = sqrt(var);
	double temp_scale;
	for(int k = 0; k< K; k++){
		if(cluster_count(k)!=0){
			temp_scale = ( nu_0 + beta.segment(P+P_two*k,P_two).squaredNorm()) /( nu_0 + P_two );
			unique_taus_sq(k) =  ((nu_0 + P_two) * temp_scale )/ rchisq_tau(rng);
			for(int p = (P+P_two*k); p<(P+P_two*(k+1));p++)
				tau_var_matrix(p) =  (1/ unique_taus_sq(k));
		}
	}

}

void FDPPSampler::store_samples(Eigen::ArrayXXd &beta_samples,
						       Eigen::ArrayXd &sigma_samples,
							   Eigen::ArrayXXd &pi_samples,
							   Eigen::ArrayXXd &tau_samples,
							   Eigen::ArrayXd &alpha_samples,
							   Eigen::ArrayXXi &cluster_assignment,
							   Eigen::ArrayXXd &yhat_samples){

	beta_samples.row(sample_ix) = beta.transpose().array();
	sigma_samples(sample_ix) = sigma;
	pi_samples.row(sample_ix) = pi.transpose();
	alpha_samples(sample_ix) = alpha;
	cluster_assignment.row(sample_ix) = iter_cluster_assignment;
	tau_samples.row(sample_ix) = sqrt(unique_taus_sq);
	yhat_samples.row(sample_ix) = yhat;
	sample_ix ++;
	int j = 0;
	for(int i = 0; i< n; i ++){
		for(int k = 0; k < K; k ++){
			if(cluster_matrix(i,k) == 1)
				P_matrix.row(i) = P_matrix.row(i) + cluster_matrix.col(k).transpose().array();
		}
	}

}

void FDPPSampler::initialize_beta(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	std::chi_squared_distribution<double> rchisq(nu_0);
	for(int k = 0; k< K; k++){
		unique_taus_sq(k) = (nu_0 / rchisq(rng));
		for(int p = (P+P_two*k); p<(P+P_two*(k+1));p++){
			beta(p) = rnorm(rng) * sigma * unique_taus_sq(k);
			tau_var_matrix(p,p) =  1/ unique_taus_sq(k);
		}
	}

}


