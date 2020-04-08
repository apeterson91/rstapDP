void GFDPPSampler::iteration_sample(std::mt19937 &rng){

	calculate_b();
	sample_cluster_labels(rng);
	update_weights(rng);

	std::uniform_real_distribution<double> runif(0,1);
	draw_z(rng);
	X_fit << Z, X_K;
	Eigen::VectorXd temp;
	diff = 1;
	int cntr = 0;
	//find mode
	Rcpp::Rcout << "cluster count: " << cluster_count << std::endl;
	Rcpp::Rcout << "initial beta: " << beta.tail(5) << std::endl;
	while(diff > 1E-6){
		temp = beta;
		eta = (X_fit * beta).array();
		W =  ((nt * exp(eta) ) / ( pow(1 + exp(eta),2) )) ; 
		mu = exp(eta) / (1 + exp(eta));
		r = eta + (y / nt - mu ) * pow(1 + exp(eta),2) / exp(eta);
		beta = ( X_fit.transpose() * W.asDiagonal() * X_fit + correction_mat ).ldlt().solve(X_fit.transpose() * W.asDiagonal() * r);
		if(cntr < 2){
			Rcpp::Rcout << " beta.tail: " << beta.tail(3) << std::endl;
			Rcpp::Rcout << " r " << (r).tail(5) << std::endl;
			Rcpp::Rcout << " W " << (W).tail(5) << std::endl;
			Rcpp::Rcout << "X * W * r " << (X_fit.transpose() * W.asDiagonal() * r).tail(5) << std::endl;
			Rcpp::Rcout << "---------------- "  << std::endl;
			//Rcpp::Rcout << "XTWX + correction: " << (X_fit.transpose() * W.asDiagonal() * X_fit + correction_mat).inverse().diagonal().tail(5) << std::endl;
		}
		diff = (beta - temp).dot(beta-temp);
		cntr ++;
		if(cntr >50){
			Rcpp::Rcout << " Uh-oh, taking a while to find mode" << std::endl;
			Rcpp::Rcout << " Current beta" << beta << std::endl;
			Rcpp::Rcout << " Current diff " << diff << std::endl;
			break;
		}
	}
	Rcpp::Rcout << "cntr: " << cntr << std::endl;
	L = ((X_fit.transpose() * W.asDiagonal() * X_fit + correction_mat).inverse()).llt().matrixL().toDenseMatrix() ;
	beta = L * z + beta;


	// error check
	if( (beta.array().isNaN() == true).any() & flag ){
		Rcpp::Rcout << "things are NaN" << std::endl;
		/*
		Rcpp::Rcout << " correction_mat " << correction_mat.diagonal() << std::endl;
		Rcpp::Rcout << "exp(eta): \n" << exp(eta).isNaN().any() << std::endl;
		Rcpp::Rcout << "(nt * exp(eta)): \n " << ( (nt * exp(eta)) ).isNaN().any() << std::endl;
		Rcpp::Rcout << "pow(1 +exp(eta),2): \n " <<  (pow(1 + exp(eta),2) ).isNaN().any() << std::endl;
		Rcpp::Rcout << "sqrt(nt * exp(eta)/...): \n " << (sqrt( (nt * exp(eta)) / pow(1 + exp(eta),2) )).isNaN().any() << std::endl;
		Rcpp::Rcout << "exp(eta): \n" << exp(eta).tail(5) << std::endl;
		Rcpp::Rcout << "(nt * exp(eta)): \n " << ( (nt * exp(eta)) ).tail(5) << std::endl;
		Rcpp::Rcout << "pow(1 +exp(eta),2): \n " <<  (pow(1 + exp(eta),2) ).tail(5) << std::endl;
		Rcpp::Rcout << "sqrt(nt * exp(eta)/...): \n " << (sqrt( (nt * exp(eta)) / pow(1 + exp(eta),2) )).tail(5)<< std::endl;
		Rcpp::Rcout << "X * beta: \n " << (X_fit * beta).tail(5) << std::endl;
		Rcpp::Rcout << "eta: \n " << (eta).tail(5) << std::endl;
		*/
		flag = false;
	}



	draw_var(rng);
}

void GFDPPSampler::draw_z(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	for(int q = 0; q < Q; q ++)
		z(q) = rnorm(rng);

}

void GFDPPSampler::draw_var(std::mt19937 &rng){

	/*
	std::chi_squared_distribution<double> rchisq_tau(nu_0 + 1);
	double temp_scale;
	for(int k = 0; k< (Q-P); k += P_two){
		for(int p_ix =0; p_ix < P_two ; p_ix ++ ){
			temp_scale = ((nu_0 + pow(beta(P + k + p_ix),2)) * W(P + k +p_ix)   ) / (1 + nu_0);
			tau_matrix(P + k+p_ix,P + k+p_ix) = sqrt(rchisq_tau(rng) / ( (nu_0+1) * temp_scale)); // inverted scale
		}
	}
	tau_var_matrix = tau_matrix * tau_matrix;
	*/

}

void GFDPPSampler::store_samples(Eigen::ArrayXXd &beta_samples,
							   Eigen::ArrayXXd &pi_samples,
							   Eigen::ArrayXXd &tau_samples,
							   Eigen::ArrayXd  &alpha_samples,
							   Eigen::ArrayXXi &cluster_assignment){

	beta_samples.row(sample_ix) = beta.transpose().array(); 
	pi_samples.row(sample_ix) = pi.transpose();
	alpha_samples(sample_ix) = alpha;
	cluster_assignment.row(sample_ix) = iter_cluster_assignment;
//	tau_samples.row(sample_ix) = Eigen::Vectord::Ones(tau_samples.cols());
	sample_ix ++;
	for(int i = 0; i< n; i ++){
		for(int j = 0; j < i ; j++)
			P_matrix(i,j) += iter_cluster_assignment(i) == iter_cluster_assignment(j) ? 1 : 0;
	}
}

void GFDPPSampler::stick_break(std::mt19937 &rng){

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

void GFDPPSampler::calculate_b(){

	b.setZero(n,K);
	int start = P;
	residual = (Z * beta.head(P)).array();
	for(int k = 0; k < K; k++){
		eta = residual + (X * beta.segment(start,P_two)).array();
		b.col(k)  = log(pi(k))  + y * eta - nt * log( 1 +  exp(eta)  );
		start += P_two;
	}
}

void GFDPPSampler::sample_cluster_labels(std::mt19937 &rng){

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
	correction_mat.setZero(Q,Q);
	for(int k = 0; k < (Q-P); k+= P_two){
		cluster_count(k_) = (iter_cluster_assignment == k_).count();
		if(cluster_count(k_) == 0)
			correction_mat.diagonal().segment(P+k,P_two) = Eigen::VectorXd::Ones(P_two);
		X_K.block(0,k,n,P_two) = cluster_matrix.col(k_).asDiagonal() * X;
		k_ ++;
	}
}

void GFDPPSampler::update_weights(std::mt19937 &rng){

	for(int k = 0; k < K; k++){
		u_posterior_beta_alpha(k) = 1 + cluster_count(k);
		u_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
	}

	stick_break(rng);
	posterior_b_alpha = 1.0 / alpha_b - (log(1 - u.array())).head(K-1).sum();
	std::gamma_distribution<double>  rgamma(posterior_a_alpha,posterior_b_alpha);
	alpha = rgamma(rng);

}
