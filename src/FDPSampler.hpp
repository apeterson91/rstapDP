#ifndef _FDPSampler_
#define _FDPSampler_

class FDPSampler
{
	private:
		const Eigen::MatrixXd &X;
		const Eigen::VectorXd &y;
		Eigen::MatrixXd S;
		Eigen::VectorXd residual;
		Eigen::MatrixXd f;
		Eigen::VectorXd beta;
		Eigen::VectorXd z;
		double sigma;
		double s;
		double var;
		int q;
		int n;
		int sample_ix = 0;

	public:
		FDPSampler(const Eigen::MatrixXd &_X,
				   const Eigen::VectorXd &_y,
				   const Eigen::VectorXd &init_beta): 
			X(_X), y(_y), beta(init_beta) 
	{
				q = X.cols(); 
				n = y.rows();
				z = Eigen::VectorXd::Zero(n);
				f = Eigen::MatrixXd::Zero(n,q);
				sigma = 1;
				var = 1;
	}

		void iteration_sample(std::mt19937 &rng);

		void calculate_residual(const int &j);

		void calculate_Smat(const int &j);

		void draw_z(std::mt19937 &rng);

		void draw_var(std::mt19937 &rng);

		void store_samples(Eigen::ArrayXXd &beta_samples,
						   Eigen::ArrayXd &sigma_samples);

		Eigen::ArrayXd get_beta() const{
			return(beta.array());
		}

		double get_sigma() const{
			return(sigma);
		}
};

#include "FDPSampler.inl"
#endif

