#include <fstream>

void FDPPSampler::iteration_sample(std::mt19937 &rng){

	calculate_b();
	sample_cluster_labels(rng);
	update_weights(rng);

	draw_z(rng);

	X_fit.block(0,P,n,Q-P) = X_K;

	V = (X_fit.transpose() * w.asDiagonal() * X_fit + PenaltyMat).inverse();

	beta = V.llt().matrixL().toDenseMatrix() * z * sigma  + 
		   V * X_fit.transpose() * w.asDiagonal() *  y; 

	if(std::isnan(beta(0)) & flag ){
		Rcpp::Rcout << "things are NaN" << std::endl;
		Rcpp::Rcout << " V block: \n" << V.block(0,0,5,5) << std::endl;
		Rcpp::Rcout << " cluster_count:\n" << cluster_count << std::endl;
		Rcpp::Rcout << " taus: \n" << unique_taus << std::endl;
		flag = false;
	}

	yhat = X_fit * beta;

	draw_var(rng);

}

void FDPPSampler::calculate_b(){

	b.setZero(n,K);
	residual = y - Z * beta.head(P);
	int start = P;
	for(int k = 0; k < K; k++){
		b.col(k)  =  log(pi(k)) - .5 * log(2 * M_PI ) - log(sigma) -  
			(.5  * precision * pow( (residual -  X  * beta.segment(start,P_two) ).array(),2));
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
	X_K.setZero(n,P_two*K);
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
	std::gamma_distribution<double> rgamma(1,1);

	for(int k = 0; k < (Q-P); k+= P_two){
		cluster_count(k_) = (iter_cluster_assignment == k_).count();
		X_K.block(0,k,n,P_two) = cluster_matrix.col(k_).asDiagonal() * X;
		if(cluster_count(k_) == 0){
			for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++ ){
				unique_taus(k_,pen_ix) = rgamma(rng);
				PenaltyMat = unique_taus(k_,pen_ix)  * S.block(0,Q*pen_ix,Q,Q);
			}
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


void FDPPSampler::draw_var(std::mt19937 &rng){

	PenaltyMat.setZero(Q,Q);
	residual = y - X_fit * beta ;
	s = (residual.transpose() * w.asDiagonal()).dot(residual) * .5;
	s += .5 * (beta.transpose() * PenaltyMat).dot(beta) ;
	std::gamma_distribution<double> rgamma(1 + n/2 + P_two, 1/(1 + s) );
	precision = rgamma(rng);
	sigma = sqrt(1 / precision);
	double temp_scale;
	for(int k = 0; k< K; k++){
		if(cluster_count(k)!=0){
			for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++){
				temp_scale = (beta.segment(P+P_two*k,P_two).transpose() * S.block(P,P+Q*pen_ix,Q-P,Q-P)).dot( 
					beta.segment(P+P_two*k,P_two)) * .5;
				std::gamma_distribution<double> rgamma_tau(1 + P_two / 2, 1/(1 + temp_scale) );
				unique_taus(k,pen_ix) = rgamma_tau(rng);
				PenaltyMat = PenaltyMat + unique_taus(k,pen_ix)  * S.block(0,Q*pen_ix,Q,Q);
			}
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
	for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++)
		tau_samples.block(sample_ix,K*pen_ix,1,K) = unique_taus.col(pen_ix);
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
	std::gamma_distribution<double> rgamma(1,1);
	for(int k = 0; k< K; k++){
		for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++ ){
			unique_taus(k,pen_ix) = rgamma(rng);
			PenaltyMat = PenaltyMat + unique_taus(k,pen_ix) * S.block(0,Q*pen_ix,Q,Q);
		}
		for(int p = (P+P_two*k); p<(P+P_two*(k+1));p++)
			beta(p) = rnorm(rng); 
	}
}


void FDPPSampler::draw_z(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	for(int q = 0; q < Q; q ++)
		z(q) = rnorm(rng);

}
