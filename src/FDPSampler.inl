void FDPSampler::iteration_sample(std::mt19937 &rng){


	calculate_f();
	for(int j = 0; j < J; j ++){
		calculate_residual(j);
		calculate_Smat(j);
		draw_z(rng);
		f_Z.col(j) =  (S * residual  + sigma * S * z);
		delta(j) = ( f_Z.col(j).array() / Z.col(j).array() )(0);
	}
	sample_cluster_labels(rng);
	update_cluster_pars(rng);

	draw_beta_var(rng);
	for(int q = 0; q < Q; q++){
		calculate_residual(q);
		sample_cluster_betas(q);
	}

	draw_var(rng);
}

void FDPSampler::calculate_f(){

	f_X = X.array() * (cluster_mat * cluster_betas.transpose()).array();
	f_Z = Z.array().rowwise() * delta.transpose().array();

}

void FDPSampler::calculate_residual(const int &j){

	residual = y;
	int p = 0;
	while( p < J){
		if(p != j)
			residual = residual - f_Z.col(p) ;
		p++;
	}
	p = 0;
	while(p < Q){
		if(p !=j)
			residual = residual - f_X.col(p);
		p++;
	}
}

void FDPSampler::calculate_Smat(const int &j){

	S = Z.col(j) * ( Z.col(j).transpose() * Z.col(j)).inverse() * Z.col(j).transpose();
}

void FDPSampler::draw_z(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	for(int i = 0; i < n; i ++)
		z(i) = rnorm(rng);
}

void FDPSampler::draw_var(std::mt19937 &rng){

	calculate_f();
	s = (y - f_X.colwise().sum() - f_Z.colwise().sum()) / ( n - q - J *K);
	std::chi_squared_distribution<double> rchisq(n - q - J * K);
	var = rchisq(rng) ;
	var = (s * n - q - J * K) / var;
	sigma = sqrt(var);
}

void FDPSampler::store_samples(Eigen::ArrayXXd &delta_samples,
						   Eigen::ArrayXd &sigma_samples){

	delta_samples.row(sample_ix) = delta.transpose().array();
	sigma_samples(sample_ix) = sigma;
	sample_ix ++;

}

void FDPSampler::stick_break(std::mt19937 &rng){

	if(initializing){
		sftrabbit::beta_distribution<> rbeta(1,alpha);
		for(int k = 0; k < (K-1) ; k++)
			u(k) = rbeta(rng);
		u(K-1) = 1;
	}else{
		for(int k = 0; k < (K-1); k++){
			sftrabbit::beta_distribution<> rbeta(u_posterior_beta_alpha(k),u_posterior_beta_beta(k));
			u(k) = rbeta(rng);
		}
		u(K-1) = 1;
	}
	for(int k = 0; k < K ; k ++)
		pi(k) = k  == 0 ? u(k) : u(k) * (1 - u.head(k)).prod();

}

void FDPSampler::calculate_b(){

	residual = y - Z*delta;
	for(int i = 0; i < n; i++){
		for(int k = 0; k < K; k++)
			b(i,k) = - .5 / var * pow(residual(i) - X.row(i).dot(cluster_betas.col(k)),2);
	}
}

void FDPSampler::sample_cluster_labels(std::mt19937 &rng){

	calculate_b();
	cluster_mat = Eigen::MatrixXd::Zeros(n,K);
	double log_factor = log(pow(10,-16)) - log(n);
	for(int i = 0; i < n; i++){
		probs = b.row(i);
		probs = probs - probs.maxCoeff();
		for(int k = 0; k < K ; k ++)
			probs(k) = probs(k) > log_factor ? exp(probs(k)) : 0 ;
		std::discrete_distribution<int> d(probs.data(), probs.data() + probs.size());
		iter_cluster_assignment(i) = d(rng);
		cluster_mat(i,iter_cluster_assignment(i)) = 1;
	}
}

void FDPSampler::update_cluster_pars(std::mt19937 &rng){

	for(int k = 0; k < K; k++)
		cluster_count(k) = (iter_cluster_assignment == k).count();

	for(int k = 0; k < K; k++){
		u_posterior_beta_alpha(k) = 1 + cluster_count(k);
		u_posterior_beta_beta(k) = alpha + cluster_count.tail(K-k-1).sum();
	}

	posterior_b_alpha = 1.0 / alpha_b - (log(1 - u.array())).head(K-1).sum();

	std::gamma_distribution<double> rgamma(posterior_a_alpha,posterior_b_alpha);
	alpha = rgamma(rng);
	stick_break(rng);

}

void FDPSampler::sample_cluster_betas(const int &q,std::mt19937 &rng){

	std::chi_squared_distribution<double> rchisq(nu_0);
	std::normal_distribution<double> rnorm(0,1);
	beta_sigma = 1 / rchisq(rng); 
	double sum;
	double mu_n;
	double var_n;
	for(int k = 0; k < K; k ++){
		if(cluster_count(k)==0)
			cluster_betas(q,k) = rnorm(rng) * beta_sigma + beta_mu;
		else{
			for(int i = 0; i < n; i ++){
				if(iter_cluster_assignment(i) == k){
					sum += residual(i);
				}
			}
			mu_n = (beta_kappa / var * beta_mu + sum / var ) / (  (beta_kappa /var) + (cluster_count(k) / var ) );
			var_n = 1 /(  ( beta_kappa / var  ) + (cluster_count(k) / var) );
			cluster_betas(q,k) = rnorm(rng) * sqrt(var_n) + mu_n;
		}
	}

}

void FDPSampler::draw_beta_var(std::mt19937 &rng){


	std::chi_squared_distribution<double> rchisq(nu_0 + n);
	beta_var = 


}
