void VI::initialize_pars(std::mt19937 &rng){

	std::normal_distribution<double> z(0,1);
	sftrabbit::beta_distribution<> rbeta(1,1); // need to correct to include alpha later
	beta.setZero(Q);
	for(int i = 0 ;i < Q; i ++)
		beta(i) = z(rng);
	tau = exp(z(rng));
	omega = exp(z(rng));
	current_iter = 0;
	current_bound = 100;
	//initialize weights
	for(int k = 0; k < (K-1); k++){
		u(k) = rbeta(rng);
	}
	u(K-1) = 1;

}

void VI::stick_break(std::mt19937 &rng){

	// for (int k = 0 ; k < K; k ++)
	//		pi(k) = k == 0 ? u(k) : u(k) * (Eigen::ArrayXd::Ones(k) - u.head(k)).prod();

}

void VI::descend_gradient(){


	V = (tau * Eigen::MatrixXd::Identity(Q,Q)  + omega * X.transpose()*X).inverse();

	// regression coefficients
	beta = omega * V * X.transpose() * y;

	omega_a_n = omega_a + N * .5;
	omega_b_n = omega_b + .5 * (y - X*beta).squaredNorm() + (X.transpose()*X * V).trace();
	
	tau_a_n = tau_a + Q * .5;
	tau_b_n = tau_b + 2 * beta.dot(beta) + V.trace();

	//residual precision
	omega = omega_a_n / omega_b_n;

	// coefficient ridge penalty
	tau = tau_a_n / tau_b_n;

}

double VI::calculate_bound(){

	double out = 0;

	out += N *.5 * boost::math::digamma(omega_a_n) - log(omega_b_n) - log(2*M_PI) ;
	out += - omega_a_n  / (2 * omega_b_n) *( y - X*beta).squaredNorm() + (X.transpose()*X*V).trace();
	out += - Q * .5 * log(2*M_PI) + Q * .5 *(boost::math::digamma(tau_a_n) - log(tau_b_n));
	out += tau_a_n / (2 * tau_b_n) * ( beta.dot(beta) + V.trace() );
	out += tau_a * log(tau_b) + (tau_a-1) *(boost::math::digamma(tau_a_n) - log(tau_b_n)) - tau_b * (tau_a_n/tau_b_n) - std::lgamma(tau_a);
	out += .5 * log(V.determinant()) + Q * .5 * (1 + log(2*M_PI));
	out += std::lgamma(tau_a_n) - (tau_a_n - 1) * boost::math::digamma(tau_a_n) - log(tau_b_n) + tau_a_n;

	return(out);
}


bool VI::has_not_converged(arrr &bounds){

	if(current_iter>max_iter){
		Rcpp::Rcout << "Max Iteration: " << current_iter << " has been reached" << std::endl;
		return(false);
	}
	double new_bound = calculate_bound();
	bounds(current_iter) = new_bound;
	current_iter ++;

	if(current_bound > new_bound){
		::Rf_error("Error, new bound should always be greater than previous lower bound ");
	}

	if(abs(current_bound - new_bound) <= 1E-6)
		return(false); // has converged, so has_not_converged = false
	else{
		current_bound = new_bound;
		return(true); // has not converged, so has_not_converged = true
	}
}

void VI::draw_pars(std::mt19937 &rng,
				  const int &num_samples,
					arrmat &beta_samples,
					arrr &sigma_samples,
					arrr &tau_samples,
					arrmat &yhat_samples){

	std::normal_distribution<double> z(0,1);
	std::gamma_distribution<double> rgamma_tau(tau_a_n,1 / tau_b_n);
	std::gamma_distribution<double> rgamma_omega(omega_a_n,1 / omega_b_n);

	for(int i = 0 ; i < num_samples; i++){
		sigma_samples(i) = sqrt( 1/ rgamma_omega(rng));
		tau_samples(i) = sqrt( 1 / rgamma_tau(rng));
		for(int j = 0; j < Q; j ++)
			beta_samples(i,j) = z(rng);
		beta_samples.row(i) = V.llt().matrixL().toDenseMatrix()* beta_samples.row(i).matrix().transpose()  + beta;
	}

}
