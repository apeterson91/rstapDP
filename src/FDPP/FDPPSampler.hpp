#ifndef _FDPPSampler_
#define _FDPPSampler_

#include "../auxiliary/beta_rng.hpp"
#include <Eigen/Dense>

class FDPPSampler
{
	private:
		const Eigen::VectorXd &y;
		const Eigen::MatrixXd &Z;
		const Eigen::MatrixXd &X;
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
		Eigen::VectorXd beta;
		Eigen::VectorXd beta_temp;
		Eigen::VectorXd z;
		Eigen::ArrayXXd b;
		Eigen::ArrayXd unique_taus;
		Eigen::ArrayXd unique_tau2s;
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
		const int n;
		const int Q;
		const int K;
		const int subset_one;
		const int subset_two;
		const int &threshold;
		int sample_ix = 0;
		int num_nonzero;
		int temp_Q;
		Eigen::MatrixXd nonzero_ics;
		bool initializing = true;
		bool flag = true;
		const bool &fix_alpha;
		const bool &logging;
		const double log_factor;

	public:
		Eigen::ArrayXXd P_matrix;
		FDPPSampler(const Eigen::VectorXd &y,
				   const Eigen::MatrixXd &Z,
				   const Eigen::MatrixXd &X,
				   const Eigen::VectorXd &w,
				   const double &alpha_a,
				   const double &alpha_b,
				   const double &tau_a,
				   const double &tau_b,
				   const double &sigma_a,
				   const double &sigma_b,
				   const int &K,
				   const int &subset_one,
				   const int &subset_two,
				   const bool &logging,
				   const int &threshold,
				   const bool &fix_alpha,
				   std::mt19937 &rng
				   ): 
			y(y), Z(Z), X(X),  w(w),
			alpha_b(alpha_b),sigma_a(sigma_a),
			sigma_b(sigma_b),
			tau_a(tau_a), tau_b(tau_b),
			P(Z.cols()), 
			P_two(X.cols()),
			n(y.rows()),
			Q(Z.cols() + X.cols()*K),
			K(K),
			subset_one(subset_one),
			subset_two(subset_two),
			threshold(threshold),
			fix_alpha(fix_alpha),
			logging(logging),
			log_factor(log(pow(10,-16)) - log(n))
	{
		num_nonzero = K;
		temp_Q = P + P_two * num_nonzero;
		PenaltyMat.setZero(Q,Q); 
		P_matrix.setZero(n,n);
		unique_taus.setZero(K);
		unique_tau2s.setZero(K);
		z.setZero(temp_Q); 
		beta.setZero(Q); 
		nonzero_ics = Eigen::MatrixXd::Identity(temp_Q,temp_Q);
		beta_temp.setZero(P + P_two*num_nonzero);
		X_K.setZero(n,P_two*num_nonzero);
		u.setZero(K);
		pi.setZero(K);
		u_posterior_beta_alpha.setZero(K);
		u_posterior_beta_beta.setZero(K);
		cluster_count.setZero(K);
		iter_cluster_assignment.setZero(n);
		probs.setZero(K);
		sigma = 1;
		precision = 1;
		std::gamma_distribution<double> rgamma(alpha_a,alpha_b);
		if(fix_alpha)
			alpha = alpha_a;
		else{
			alpha = rgamma(rng);
			posterior_a_alpha = alpha_a + K - 1;
		}
		stick_break(rng);
		cluster_matrix.setZero(n,K);
		initialize_beta(rng);
		initializing = false;
	}

		void iteration_sample(std::mt19937 &rng);

		void draw_z(std::mt19937 &rng);

		void draw_var(std::mt19937 &rng);

		void store_samples(Eigen::ArrayXXd &beta_samples,
						   Eigen::ArrayXd &sigma_samples,
						   Eigen::ArrayXXd &pi_samples,
						   Eigen::ArrayXXd &tau_samples,
						   Eigen::ArrayXXd &tau2_samples,
						   Eigen::ArrayXd &alpha_samples,
						   Eigen::ArrayXXi &cluster_assignment,
						   Eigen::ArrayXXd &yhat_samples);

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

		void update_penaltymat(const int &k, const int &subset, const bool second);

		void adjust_zero_clusters(std::mt19937 &rng);

		double calculate_penalty_scale(const int &k,const int &subset, const bool second);

		void adjust_beta(std::mt19937 &rng);

		void log_message(const std::string& input){
			if(logging)
				Rcpp::Rcout << "Log: " <<  input << std::endl;
		}

};

#include "FDPPSampler.inl"
#endif

