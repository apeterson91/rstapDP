#ifndef _FDPSampler_
#define _FDPSampler_

class FDPSampler
{
	private:
		const Eigen::MatrixXd &X;
		const Eigen::VectorXd &y;
		const Eigen::MatrixXd &Z;
		Eigen::MatrixXd S;
		Eigen::VectorXd residual;
		Eigen::MatrixXd f_X;
		Eigen::MatrixXd f_Z;
		Eigen::VectorXd delta;
		Eigen::MatrixXd cluster_betas;
		Eigen::VectorXd z;
		Eigen::ArrayXi iter_cluster_assignment;
		Eigen::MatrixXd cluster_mat;
		Eigen::ArrayXXd b;
		Eigen::ArrayXd pi;
		Eigen::ArrayXd probs;
		Eigen::ArrayXd u;
		Eigen::ArrayXi dp_ics;
		Eigen::ArrayXd u_posterior_beta_alpha;
		Eigen::ArrayXd u_posterior_beta_beta; 
		Eigen::ArrayXd cluster_count(k);
		double sigma;
		double s;
		double var;
		double posterior_a_alpha;
		double posterior_b_alpha;
		double alpha;
		const double &b_alpha;
		const double &beta_mu;
		const double &beta_kappa;
		double beta_sigma;
		int q;
		int n;
		int J;
		int sample_ix = 0;
		const int &K;
		bool initializing = true;

	public:
		FDPSampler(const Eigen::MatrixXd &_X,
				   const Eigen::MatrixXd &_Z,
				   const Eigen::VectorXd &_y,
				   const int &_K,
				   const Eigen::VectorXd &init_delta,
				   const Eigen::ArrayXi &_dp_ics,
				   const double &a_alpha,
				   const double &b_alpha,
				   const double &beta_kappa,
				   const double &beta_mu,
				   std::mt19937 &rng): 
			X(_X), y(_y), delta(init_delta), dp_ics(_dp_ics), 
			K(_K), b_alpha(b_alpha), beta_mu(beta_mu),
			beta_kappa(beta_kappa)
	{
		Q = X.cols(); 
		J = Z.cols();
		n = y.rows();
		posterior_a_alpha =  a_alpha + K - 1;
		z = Eigen::VectorXd::Zero(n);
		f_X = Eigen::MatrixXd::Zero(n,Q);
		f_Z = Eigen::MatrixXd::Zero(n,J);
		delta = Eigen::VectorXd::Zero(J);
		cluster_betas = Eigen::MatrixXd(Q,K);
		b = Eigen::ArrayXXd::Zero(n,K);
		sigma = 1;
		var = 1;
		stick_break(rng);
		sample_cluster_labels(rng);
		initializing = false;
	}

		void iteration_sample(std::mt19937 &rng);

		void calculate_f();

		void calculate_residual(const int &j);

		void calculate_Smat(const int &j);

		void draw_z(std::mt19937 &rng);

		void draw_var(std::mt19937 &rng);

		void store_samples(Eigen::ArrayXXd &delta_samples,
						   Eigen::ArrayXd &sigma_samples);

		void stick_break(std::mt19937 &rng);

		void sample_cluster_betas(const int &q);

		Eigen::ArrayXd get_delta() const{
			return(delta.array());
		}

		double get_sigma() const{
			return(sigma);
		}
};

#include "FDPSampler.inl"
#endif

