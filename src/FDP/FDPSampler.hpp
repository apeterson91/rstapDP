#ifndef _FDPSampler_
#define _FDPSampler_

#include "../auxiliary/beta_rng.hpp"

class FDPSampler
{
	private:
		const Eigen::MatrixXd &X;
		const Eigen::MatrixXd &Z;
		const Eigen::VectorXd &y;
		const Eigen::VectorXd &w;
		Eigen::MatrixXd X_fit;
		Eigen::MatrixXd X_K;
		Eigen::MatrixXd V;
		Eigen::MatrixXd correction_mat;
		Eigen::ArrayXd pi;
		Eigen::ArrayXd u;
		Eigen::ArrayXi iter_cluster_assignment;
		Eigen::MatrixXd cluster_matrix;
		Eigen::ArrayXi cluster_count;
		Eigen::ArrayXd probs;
		Eigen::VectorXd residual;
		Eigen::VectorXd beta;
		Eigen::VectorXd z;
		Eigen::ArrayXXd b;
		Eigen::MatrixXd tau_matrix;
		const double alpha_b;
		const Eigen::ArrayXd tau_0;
		Eigen::ArrayXd u_posterior_beta_alpha;
		Eigen::ArrayXd u_posterior_beta_beta;
		double sigma;
		double s;
		double var;
		double alpha;
		double posterior_a_alpha;
		double posterior_b_alpha;
		const int P;
		const int P_two;
		const int n;
		const int K;
		const int Q;
		int sample_ix = 0;
		bool initializing = true;

	public:
		Eigen::MatrixXd P_matrix;
		FDPSampler(const Eigen::VectorXd &y,
				   const Eigen::MatrixXd &Z,
				   const Eigen::MatrixXd &X,
				   const Eigen::VectorXd &w,
				   const double &alpha_a,
				   const double &alpha_b,
				   const Eigen::ArrayXd &tau_0,
				   const int &K,
				   std::mt19937 &rng
				   ): 
			X(X), Z(Z), y(y),w(w),
			alpha_b(alpha_b), tau_0(tau_0),
			n(y.rows()), P(Z.cols()), K(K),
			Q(Z.cols() + X.cols()*K), P_two(X.cols())
	{

		P_matrix.setZero(n,n);
		correction_mat.setZero(Q,Q);
		tau_matrix.setZero(Q,Q);
		X_fit.setZero(n,Q);
		z.setZero(Q); 
		beta.setZero(Q); 
		X_K.setZero(n,P_two*K);
		u.setZero(K);
		pi.setZero(K);
		u_posterior_beta_alpha.setZero(K);
		u_posterior_beta_beta.setZero(K);
		cluster_count.setZero(K);
		iter_cluster_assignment.setZero(n);
		probs.setZero(K);
		sigma = 1;
		var = 1;
		std::gamma_distribution<double> rgamma(alpha_a,alpha_b);
		alpha = rgamma(rng);
		posterior_a_alpha = alpha_a + K - 1;
		stick_break(rng);
		cluster_matrix.setZero(n,K);
		initializing = false;
		initialize_beta(rng);

	}

		void iteration_sample(std::mt19937 &rng);

		void calculate_residual(const int &j);

		void calculate_Smat(const int &j);

		void draw_z(std::mt19937 &rng);

		void draw_var(std::mt19937 &rng);

		void store_samples(Eigen::ArrayXXd &beta_samples,
						   Eigen::ArrayXd &sigma_samples,
						   Eigen::ArrayXXd &tau_samples,
						   Eigen::ArrayXXd &pi_samples,
						   Eigen::ArrayXd &alpha_samples,
						   Eigen::ArrayXXi &cluster_assignment);

		Eigen::ArrayXd get_beta() const{
			return(beta.array());
		}

		double get_sigma() const{
			return(sigma);
		}

		void stick_break(std::mt19937 &rng);

		void calculate_b();

		void sample_cluster_labels(std::mt19937 &rng);

		void update_weights(std::mt19937 &rng);

		void initialize_beta(std::mt19937 &rng);

};

#include "FDPSampler.inl"
#endif

