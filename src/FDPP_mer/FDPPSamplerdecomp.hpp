#ifndef _FDPPMERD_
#define _FDPPMERD_

#include "../auxiliary/beta_rng.hpp"
#include "../auxiliary/wishart_rng.hpp"
#include <Eigen/Dense>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::ArrayXXd arr;

class FDPPSamplerdecomp 
{
	private:
		const Eigen::VectorXd &y;
		const Eigen::MatrixXd &Z;
		const Eigen::MatrixXd &X_b;
		const Eigen::MatrixXd &X_w;
		const arr W;
		const SpMat &subj_mat;
		const Eigen::MatrixXd &S_b;
		const Eigen::MatrixXd &S_w;
		const Eigen::VectorXd &w;
		Eigen::VectorXd yhat;
		Eigen::MatrixXd X_fit;
		Eigen::MatrixXd X_K;
		Eigen::MatrixXd V;
		Eigen::MatrixXd PenaltyMat;
		Eigen::ArrayXd pi;
		Eigen::ArrayXd u;
		Eigen::ArrayXi iter_cluster_assignment;
		Eigen::MatrixXd cluster_matrix;
		Eigen::ArrayXi cluster_count;
		Eigen::ArrayXd probs;
		Eigen::VectorXd residual;
		Eigen::VectorXd Wb;
		Eigen::VectorXd beta;
		Eigen::VectorXd beta_temp;
		Eigen::VectorXd z;
		Eigen::VectorXd z_b;
		Eigen::ArrayXXd b;
		Eigen::MatrixXd subj_b;
		Eigen::MatrixXd subj_D;
		const int subj_D_df;
		const Eigen::ArrayXi subj_n;
		Eigen::ArrayXXd unique_taus_b;
		Eigen::ArrayXXd unique_taus_w;
		const double alpha_b;
		const double sigma_a; 
		const double sigma_b; 
		const double tau_a; 
		const double tau_b;
		Eigen::ArrayXd u_posterior_beta_alpha;
		Eigen::ArrayXd u_posterior_beta_beta;
		double sigma;
		double s;
		double precision;
		double alpha;
		double posterior_a_alpha;
		double posterior_b_alpha;
		const int P;
		const int P_two;
		const int P_b;
		const int P_w;
		const int N;
		const int n;
		const int K;
		const int Q;
		const int num_penalties;
		int sample_ix = 0;
		int num_nonzero;
		int temp_Q;
		Eigen::MatrixXd nonzero_ics;
		bool initializing = true;
		bool flag = true;
		const bool fix_alpha;
		const double log_factor;
	public:
		Eigen::ArrayXXd P_matrix;
		FDPPSamplerdecomp(
				const Eigen::VectorXd &y,
				const Eigen::MatrixXd &Z,
				const Eigen::MatrixXd &X_b,
				const arr &W,
				const Eigen::MatrixXd &S_b,
				const Eigen::VectorXd &w,
				const SpMat &subj_mat,
				const Eigen::ArrayXi &subj_n,
				const double &alpha_a,
				const double &alpha_b,
				const double &sigma_a,
				const double &sigma_b,
				const double &tau_a,
				const double &tau_b,
				const int &K,
				const int &num_penalties,
				const bool &fix_alpha,
				std::mt19937 &rng,
				const Eigen::MatrixXd &X_w,
				const Eigen::MatrixXd &S_w) : 
			y(y), Z(Z), X_b(X_b), X_w(X_w),
			W(W), subj_mat(subj_mat),
			S_b(S_b),  S_w(S_w),w(w),
			subj_D_df(subj_mat.cols()-W.cols()+1),
			subj_n(subj_n),
			alpha_b(alpha_b),sigma_a(sigma_a),
			sigma_b(sigma_b),
			tau_a(tau_a), tau_b(tau_b),
			P(Z.cols()), 
			P_two(X_b.cols() + X_w.cols()),
			P_b(X_b.cols()),
			P_w(X_w.cols()),
			N(y.rows()), n(subj_mat.cols()),
			K(K),
			Q(Z.cols() + (X_w.cols() +  X_b.cols())*K), 
			num_penalties(num_penalties),
			fix_alpha(fix_alpha),
			log_factor(log(pow(10,-16)) - log(N))
	{
		std::gamma_distribution<double> rgamma(alpha_a,alpha_b);
		if(fix_alpha)
			alpha = alpha_a;
		else{
			alpha = rgamma(rng);
			posterior_a_alpha = alpha_a + K - 1;
		}
		initialization(rng);
	}

		void iteration_sample(std::mt19937 &rng);

		void draw_z(std::mt19937 &rng);

		void draw_zb(std::mt19937 &rng);

		void draw_var(std::mt19937 &rng);

		void store_samples(Eigen::ArrayXXd &beta_samples,
						   Eigen::ArrayXd &sigma_samples,
						   Eigen::ArrayXXd &pi_samples,
						   Eigen::ArrayXXd &tau_samples_b,
						   Eigen::ArrayXXd &tau_samples_w,
						   Eigen::ArrayXd &alpha_samples,
						   Eigen::ArrayXXi &cluster_assignment,
						   Eigen::ArrayXXd &yhat_samples,
						   Eigen::ArrayXXd &subj_b_samples,
						   Eigen::ArrayXXd &subj_D_samples
						   );

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

		void update_penaltymat(const int &k, const int &pen_ix);

		void adjust_zero_clusters(std::mt19937 &rng);

		double calculate_penalty_scale(const int &k, const int &pen_ix, const bool &between);

		void adjust_beta(std::mt19937 &rng);

		void draw_subj_b(std::mt19937 &rng);

		void draw_subj_D(std::mt19937 &rng);

		void calculate_Wb();

		void initialization(std::mt19937 &rng);


};

#include "FDPPSamplerdecomp.inl"
#endif
