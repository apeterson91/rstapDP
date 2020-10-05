// inline functions for FDPPSamplerdecomp class

void FDPPSamplerdecomp::iteration_sample(std::mt19937 &rng){


	// update cluster labels
	calculate_b();
	sample_cluster_labels(rng);
	update_weights(rng);

	// create appropriately sized containers
	beta_temp.setZero(temp_Q);
	X_fit.setZero(N,temp_Q);
	X_fit.block(0,0,N,P) = Z;
	draw_z(rng);
	X_fit.block(0,P,N,P_two*num_nonzero) = X_K;

	// calculate solution
	V = (X_fit.transpose() * w.asDiagonal() * X_fit + 
			nonzero_ics.transpose() * PenaltyMat * nonzero_ics ).inverse();

	beta_temp = V.llt().matrixL().toDenseMatrix() * sigma * z + 
		V * X_fit.transpose() * w.asDiagonal() * (y-Wb) ;

	beta = nonzero_ics * beta_temp;

	adjust_beta(rng);

	// check for errors
	if(std::isnan(beta(0)) & flag ){
		Rcpp::Rcout << "things are NaN" << std::endl;
		Rcpp::Rcout << " V block: \n" << V.block(0,0,5,5) << std::endl;
		flag = false;
	}

	// update subject terms
	draw_subj_D(rng);

	draw_subj_b(rng);

	// draw variance terms
	draw_var(rng);

	calculate_Wb();

	// calculate predicted value
	yhat = X_fit * beta + Wb;
}

void FDPPSamplerdecomp::calculate_Wb(){
	
	Wb = ((W * (subj_mat * subj_b).array())).rowwise().sum().matrix();
}

void FDPPSamplerdecomp::calculate_b(){

	b.setZero(n,K);
	residual = y - Z * beta.head(P) - Wb;
	int start_b = P;
	int start_w = P + P_b;
	Eigen::VectorXd temp;
	for(int k = 0; k < K; k++){
		temp = ( .5 * log(2 * M_PI ) - log(sigma * (1/w.array()))  -  
			(.5  *  precision * w.array() * pow( (residual -  X_b  * beta.segment(start_b,P_b) - X_w * beta.segment(start_w,P_w) ).array(),2))).matrix(); 
		// sum across measurements 
		b.col(k)  =  log(pi(k)) + (subj_mat.transpose() * temp   ).array();
		start_b += P_two;
		start_w += P_two;
	}

}

void FDPPSamplerdecomp::stick_break(std::mt19937 &rng){

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

void FDPPSamplerdecomp::sample_cluster_labels(std::mt19937 &rng){

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

	adjust_zero_clusters(rng);

}

void FDPPSamplerdecomp::adjust_zero_clusters(std::mt19937 &rng){

	std::gamma_distribution<double> rgamma(1,1);
	num_nonzero = K;

	for(int k = 0; k < K; k++){
		cluster_count(k) = (iter_cluster_assignment == k).count();
		if(cluster_count(k) == 0)
			num_nonzero --;
	}

	if(num_nonzero<=0)
		throw std::range_error("cannot have <= 0 nonzero clusters - error present");

	temp_Q = P + P_two * num_nonzero;
	nonzero_ics.setZero(Q,temp_Q);
	nonzero_ics.block(0,0,P,P) = Eigen::MatrixXd::Identity(P,P);
	X_K.setZero(N,P_two*num_nonzero);

	int k_ = 0;
	Eigen::MatrixXd temp;
	for(int k = 0; k < K; k++){
		if(cluster_count(k)!=0){
			temp = (subj_mat * cluster_matrix.col(k).asDiagonal() * subj_mat.transpose() );
			X_K.block(0,k_*P_two,N,P_b) = temp.diagonal().asDiagonal() * X_b;
			X_K.block(0,P_b + k_*P_two,N,P_w) = temp.diagonal().asDiagonal() * X_w;
			nonzero_ics.block(P+k*P_two,P+k_*P_two,P_two,P_two) = Eigen::MatrixXd::Identity(P_two,P_two);
			k_ ++;
		}
	}
}

void FDPPSamplerdecomp::update_weights(std::mt19937 &rng){

	stick_break(rng);
	if(!fix_alpha){
		posterior_b_alpha = 1.0 / alpha_b - (log(1 - u.array())).head(K-1).sum();
		posterior_b_alpha = 1.0 / posterior_b_alpha;
		std::gamma_distribution<double>  rgamma(posterior_a_alpha,posterior_b_alpha);
		alpha = rgamma(rng);
	}
}


void FDPPSamplerdecomp::draw_var(std::mt19937 &rng){

	residual = y - X_fit * beta_temp - Wb;
	s = (residual.transpose() * w.asDiagonal()).dot(residual) * .5;
	s += .5 * (beta_temp.transpose() * nonzero_ics.transpose() * PenaltyMat * nonzero_ics).dot(beta_temp) ;
	std::gamma_distribution<double> rgamma(sigma_a + N/2 + P_two, 1/( (1 / sigma_b) + s) );
	precision = rgamma(rng);
	sigma = sqrt(1 / precision);
	double temp_scale;
	PenaltyMat.setZero(Q,Q);
	for(int k = 0; k< K; k++){
		if(cluster_count(k)!=0){
			for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++){
				temp_scale = calculate_penalty_scale(k,pen_ix,true);
				std::gamma_distribution<double> rgamma_tau_b(tau_a + P_two / 2, 1/( (1/tau_b) + temp_scale) );
				unique_taus_b(k,pen_ix) = rgamma_tau_b(rng);
				// repeat for unique_taus_w
				temp_scale = calculate_penalty_scale(k,pen_ix,false);
				std::gamma_distribution<double> rgamma_tau_w(tau_a + P_two / 2, 1/( (1/tau_b) + temp_scale) );
				unique_taus_w(k,pen_ix) = rgamma_tau_w(rng);
				update_penaltymat(k,pen_ix);
			}
		}else{
			std::gamma_distribution<double> rgamma_tau_prior(tau_a,1/tau_b);
			for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++){
				unique_taus_b(k,pen_ix) = rgamma_tau_prior(rng);
				unique_taus_w(k,pen_ix) = rgamma_tau_prior(rng);
				update_penaltymat(k,pen_ix);
			}
		}
	}

}

double FDPPSamplerdecomp::calculate_penalty_scale(const int &k, const int &pen_ix,const bool &between){

	double out = 0;
	int beta_ix = P;
	int S_rix = P;
	int S_cix = P; 
	if(between){
		beta_ix += P_two *k;
		S_rix += P_b*k;
		S_cix += pen_ix * (P_b *K + P) + P_b*k;
		out = (beta.segment(beta_ix,P_b).transpose() * S_b.block(S_rix,S_cix,P_b,P_b)).dot(
				beta.segment(beta_ix,P_b));
	}else{
		beta_ix += P_b + P_two*k;
		S_rix += P_w*k;
		S_cix += pen_ix * (P_w*K + P) +  P_w*k;
		out = (beta.segment(beta_ix,P_w).transpose() * S_w.block(S_rix,S_cix,P_w,P_w)).dot(
				beta.segment(beta_ix,P_w));
	}
	out *= .5;
	return(out);
}

void FDPPSamplerdecomp::store_samples(Eigen::ArrayXXd &beta_samples,
								   Eigen::ArrayXd &sigma_samples,
								   Eigen::ArrayXXd &pi_samples,
								   Eigen::ArrayXXd &tau_samples_b,
								   Eigen::ArrayXXd &tau_samples_w,
								   Eigen::ArrayXd &alpha_samples,
								   Eigen::ArrayXXi &cluster_assignment,
								   Eigen::ArrayXXd &yhat_samples,
								   Eigen::ArrayXXd &subj_b_samples,
								   Eigen::ArrayXXd &subj_D_samples
							   ){

	beta_samples.row(sample_ix) = beta.transpose().array();
	sigma_samples(sample_ix) = sigma;
	pi_samples.row(sample_ix) = pi.transpose();
	alpha_samples(sample_ix) = alpha;
	cluster_assignment.row(sample_ix) = iter_cluster_assignment;
	for(int i = 0 ; i < W.cols() ; i ++)
		subj_b_samples.block(sample_ix,n*i,1,n) = subj_b.col(i).transpose();
	Eigen::MatrixXd temp = subj_D.inverse();
	Eigen::Map<Eigen::RowVectorXd> subj_D_row(temp.data(),temp.size());
	subj_D_samples.row(sample_ix) = subj_D_row;
	for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++){
		tau_samples_b.block(sample_ix,K*pen_ix,1,K) = unique_taus_b.col(pen_ix);
		tau_samples_w.block(sample_ix,K*pen_ix,1,K) = unique_taus_w.col(pen_ix);
	}
	yhat_samples.row(sample_ix) = yhat;
	sample_ix ++;
	for(int i = 0; i< n; i ++){
		for(int k = 0; k < K; k ++){
			if(cluster_matrix(i,k) == 1)
				P_matrix.row(i) = P_matrix.row(i) + cluster_matrix.col(k).transpose().array();
		}
	}

}

void FDPPSamplerdecomp::initialize_beta(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);
	std::gamma_distribution<double> rgamma(tau_a,(1/tau_b));
	for(int p=0;p<Q;p++)
		beta(p) = rnorm(rng);

	for(int k = 0; k< K; k++){
		for(int pen_ix = 0; pen_ix < num_penalties; pen_ix ++ ){
			unique_taus_b(k,pen_ix) = rgamma(rng);
			unique_taus_w(k,pen_ix) = rgamma(rng);
			update_penaltymat(k,pen_ix);
		}
	}

	for(int i = 0; i < n ; i ++){
		for(int q = 0; q < W.cols() ; q ++)
			subj_b(i,q) = rnorm(rng);
	}

}

void FDPPSamplerdecomp::draw_z(std::mt19937 &rng){

	z.setZero(temp_Q);
	std::normal_distribution<double> rnorm(0,1);
	for(int q = 0; q < temp_Q; q ++)
		z(q) = rnorm(rng);

}

void FDPPSamplerdecomp::draw_zb(std::mt19937 &rng){

	z_b.setZero(W.cols());
	std::normal_distribution<double> rnorm(0,1);
	for(int ix = 0; ix < W.cols(); ix ++)
		z_b(ix) = rnorm(rng);
}

void FDPPSamplerdecomp::update_penaltymat(const int &k,const int &pen_ix){

	int row_ix = P + P_two*k;
	int Srow_ix = P + P_b*k;
	int Scol_ix = P + P_b*k + (P_b*K + P)*pen_ix;
	// between
	PenaltyMat.block(row_ix,row_ix,P_b,P_b) = 
				PenaltyMat.block(row_ix,row_ix,P_b,P_b) + 
				unique_taus_b(k,pen_ix) * S_b.block(Srow_ix,Scol_ix,P_b,P_b);
	// within
	row_ix = P + P_two*k + P_b;
	Srow_ix = P + P_w*k;
	Scol_ix = P + P_w*k + (P_w*K + P)*pen_ix;
	PenaltyMat.block(row_ix,row_ix,P_w,P_w) = 
				PenaltyMat.block(row_ix,row_ix,P_w,P_w) + 
				unique_taus_w(k,pen_ix) * S_w.block(Srow_ix,Scol_ix,P_w,P_w);
}

void FDPPSamplerdecomp::adjust_beta(std::mt19937 &rng){

	std::normal_distribution<double> rnorm(0,1);

	for(int k = 0; k < K ; k++){
		if(cluster_count(k)==0){
			for(int p = 0; p < P_b; p++)
				beta(P+P_two*k+p) = rnorm(rng)*sigma * sqrt(1/unique_taus_b(k,0));
			for(int p = 0; p < P_w; p++)
				beta(P+P_two*k+P_b+p) = rnorm(rng)*sigma * sqrt(1/unique_taus_w(k,0));
		}
	}
}

void FDPPSamplerdecomp::draw_subj_b(std::mt19937 &rng){

	residual = y - X_fit * beta_temp;
	int row_ix = 0;
	int num_cols = W.cols();
	Eigen::MatrixXd temp;
	Eigen::VectorXd temp_res;
	for(int i =0 ; i < n ; i++){
		draw_zb(rng);
		temp = W.block(row_ix,0,subj_n(i),num_cols);
		temp_res = residual.segment(row_ix,subj_n(i));
		subj_b.row(i) = ((temp.transpose() * temp + sigma*subj_D).inverse() * temp.transpose() * temp_res + 
			(temp.transpose() * temp + sigma * subj_D).inverse().llt().matrixL().toDenseMatrix() * z_b).transpose();
		if((i+1)< n)
			row_ix += subj_n(i);
	}

}

void FDPPSamplerdecomp::draw_subj_D(std::mt19937 &rng){

	Eigen::MatrixXd V;

	V = subj_b.transpose() * subj_b;

	subj_D = draw_wishart(V,subj_D_df,rng);

}

void FDPPSamplerdecomp::initialization(std::mt19937 &rng){

	num_nonzero = K;
	temp_Q = P + P_two * num_nonzero;
	PenaltyMat.setZero(Q,Q); 
	P_matrix.setZero(n,n);
	unique_taus_b.setZero(K,num_penalties);
	unique_taus_w.setZero(K,num_penalties);
	b.setZero(n,K);
	z.setZero(temp_Q); 
	z_b.setZero(W.cols()); 
	beta.setZero(Q); 
	subj_D = Eigen::MatrixXd::Identity(W.cols(),W.cols());
	subj_b.setZero(n,W.cols());
	nonzero_ics = Eigen::MatrixXd::Identity(temp_Q,temp_Q);
	beta_temp.setZero(P + P_two*num_nonzero);
	X_K.setZero(N,P_two*num_nonzero);
	u.setZero(K);
	pi.setZero(K);
	u_posterior_beta_alpha.setZero(K);
	u_posterior_beta_beta.setZero(K);
	cluster_count.setZero(K);
	iter_cluster_assignment.setZero(n);
	probs.setZero(K);
	sigma = 1;
	precision = 1;
	stick_break(rng);
	cluster_matrix.setZero(n,K);
	initialize_beta(rng);
	initializing = false;
	if(y.rows() != N)
		throw std::range_error("length(y) != N");

	if(Z.rows() != N)
		throw std::range_error("nrow(Z) != N");
	if(X_b.rows() != N)
		throw std::range_error("nrow(X) != N");
	if(X_w.rows() != N)
		throw std::range_error("nrow(X) != N");

	if(subj_mat.rows() != N)
		throw std::range_error("nrow(subj_mat) != N");

	if(subj_mat.cols() != n)
		throw std::range_error("ncol(subj_mat) != n");
	
	if(b.rows() != n ){
		Rcpp::Rcout << "b rows: " << b.rows() << std::endl;
		Rcpp::Rcout << "n     : " << n << std::endl;
		throw std::range_error("nrow(b) != n");
	}

	if(subj_mat.cols() != subj_b.rows()){
		Rcpp::Rcout << "subj_mat.cols(): " << subj_mat.cols() << std::endl;
		Rcpp::Rcout << "subj_b.rows(): " << subj_b.rows() << std::endl;
		throw std::range_error("nrow(b) != n");
	}
	calculate_Wb();

}
