void FDPSampler::iteration_sample(std::mt19937 &rng){

	f = X.array().rowwise() * beta.transpose().array();
	for(int j = 0; j < q; j ++){
		calculate_residual(j);
		calculate_Smat(j);
		draw_z(rng);
		f.col(j) =  (S * residual  + sigma * S * z);
		beta(j) = (f.col(j).array() /X.col(j).array())(0);
	}

	draw_var(rng);
}

void FDPSampler::calculate_residual(const int &j){

	residual = y;
	int p = 0;
	while( p < q){
		if(p != j)
			residual = residual - f.col(p) ;
		p++ ;
	}
}

void FDPSampler::calculate_Smat(const int &j){

	S = X.col(j) * ( X.col(j).transpose() * X.col(j)).inverse() * X.col(j).transpose();
}

void FDPSampler::draw_z(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	for(int i = 0; i < n; i ++)
		z(i) = rnorm(rng);
}

void FDPSampler::draw_var(std::mt19937 &rng){

	s = (y - X*beta).dot(y-X*beta) / (n - q);
	std::chi_squared_distribution<double> rchisq(n - q);
	var = rchisq(rng) ;
	var = (s * n-q) / var;
	sigma = sqrt(var);
}

void FDPSampler::store_samples(Eigen::ArrayXXd &beta_samples,
						   Eigen::ArrayXd &sigma_samples){

	beta_samples.row(sample_ix) = beta.transpose().array();
	sigma_samples(sample_ix) = sigma;
	sample_ix ++;

}
