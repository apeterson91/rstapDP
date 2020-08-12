#ifndef _WISHARTRNG_
#define _WISHARTRNG_


Eigen::MatrixXd draw_wishart(const Eigen::MatrixXd &V,const int &df,std::mt19937 &rng){
	
	Eigen::MatrixXd L;
	L = V.llt().matrixL().toDenseMatrix();

	Eigen::MatrixXd A;
	A.setZero(V.rows(),V.cols());

	std::normal_distribution<double> rnorm(0,1);

	for(int row_ix = 0; row_ix < A.rows(); row_ix++){
		for(int col_ix = 0; col_ix <= row_ix; col_ix++){
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
