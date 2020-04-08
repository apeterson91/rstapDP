#ifndef _GFDPPSampler_
#define _GFDPPSampler_

#include "beta_rng.hpp"

class GFDPPSampler
{
	private:
		const Eigen::MatrixXd &X;
		const Eigen::MatrixXd &Z;
		const Eigen::ArrayXd &y;
		const Eigen::ArrayXd &nt;
		Eigen::VectorXd beta;
		Eigen::VectorXd temp;
		Eigen::VectorXd r;
		Eigen::ArrayXd residual;
		Eigen::ArrayXd eta;
		Eigen::ArrayXd mu;
		Eigen::VectorXd W;
		Eigen::MatrixXd L;
		Eigen::MatrixXd X_fit;
		Eigen::MatrixXd X_K;
		Eigen::ArrayXd pi;
		Eigen::ArrayXd u;
		Eigen::ArrayXi iter_cluster_assignment;
		Eigen::MatrixXd cluster_matrix;
		Eigen::ArrayXi cluster_count;
		Eigen::ArrayXd probs;
		Eigen::VectorXd z;
		Eigen::ArrayXXd b;
		Eigen::MatrixXd correction_mat ;
		Eigen::ArrayXd u_posterior_beta_alpha;
		Eigen::ArrayXd u_posterior_beta_beta;
		Eigen::ArrayXd temp_f;
		double alpha;
		double posterior_a_alpha;
		double posterior_b_alpha;
		const double alpha_b;
		const double nu_0;
		double sigma;
		const int P;
		const int P_two;
		const int n;
		const int K;
		const int Q;
		int sample_ix = 0;
		bool initializing = true;
		bool flag = true;
		double diff;
	public:
		Eigen::MatrixXd P_matrix;
		GFDPPSampler(
				   const Eigen::ArrayXd &y,
				   const Eigen::ArrayXd &nt,
				   const Eigen::MatrixXd &Z,
				   const Eigen::MatrixXd &X,
				   const double &alpha_a,
				   const double &alpha_b,
				   const double &nu_0,
				   const int &K,
				   std::mt19937 &rng
				   ): 
			X(X), Z(Z), y(y), nt(nt),
			alpha_b(alpha_b),nu_0(nu_0),
			n(y.rows()), P(Z.cols()), K(K),
			Q(Z.cols() + X.cols()*K), P_two(X.cols())
	{

		//set a bunch of things to zero
		eta.setZero(n);
		P_matrix.setZero(n,n);
		X_fit.setZero(n,Q);
		z.setZero(n); 
		X_K.setZero(n,P_two*K);
		u.setZero(K);
		pi.setZero(K);
		beta.setZero(Q);
		u_posterior_beta_alpha.setZero(K);
		u_posterior_beta_beta.setZero(K);
		cluster_count.setZero(K);
		iter_cluster_assignment.setZero(n);
		cluster_matrix.setZero(n,K);
		probs.setZero(K);
	correction_mat.setZero(Q,Q);
		//initialize  vars
		std::gamma_distribution<double> rgamma(alpha_a,alpha_b);
		alpha = rgamma(rng);
		posterior_a_alpha = alpha_a + K - 1;
		stick_break(rng);
		W.diagonal() = Eigen::VectorXd::Ones(n);
		initializing = false;
	}

		void iteration_sample(std::mt19937 &rng);


		void draw_z(std::mt19937 &rng);

		void draw_var(std::mt19937 &rng);

		void store_samples(Eigen::ArrayXXd &beta_samples,
						   Eigen::ArrayXXd &pi_samples,
						   Eigen::ArrayXXd &tau_samples,
						   Eigen::ArrayXd &alpha_samples,
						   Eigen::ArrayXXi &cluster_assignment);

		void stick_break(std::mt19937 &rng);

		void calculate_b();

		void sample_cluster_labels(std::mt19937 &rng);

		void update_weights(std::mt19937 &rng);


};

#include "GFDPPSampler.inl"
#endif
