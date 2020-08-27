#ifndef _WISHARTRNG_
#define _WISHARTRNG_

// Generated Random Wishart Matrices
// 
// Draws a random matrix from the Wishart distribution with df degrees of freedom and
// scale matrix V.
// @param V scale matrix
// @param df degrees of freedom for wishart distribution
// @param rng random number generator engine
// @return Eigen::MatrixXd Wishart draw
Eigen::MatrixXd draw_wishart(const Eigen::MatrixXd &V,const int &df,std::mt19937 &rng){
	
	Eigen::MatrixXd L;
	L = V.inverse().llt().matrixL().toDenseMatrix();

	Eigen::MatrixXd A;
	A.setZero(V.rows(),V.cols());

	std::normal_distribution<double> rnorm(0,1);

	for(unsigned int row_ix = 0; row_ix < A.rows(); row_ix++){
		for(unsigned int col_ix = 0; col_ix <= row_ix; col_ix++){
			if(col_ix==row_ix){
				std::chi_squared_distribution<double> rchisq(df-row_ix);
				A(row_ix,col_ix) = sqrt(rchisq(rng));
			}else{
				A(row_ix,col_ix) = rnorm(rng);
			}
		}
	}

	
	Eigen::MatrixXd out(V.rows(),V.cols());

	out = L * A * A.transpose() * L.transpose();

	return(out);
}


#endif
