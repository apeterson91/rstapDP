#ifndef _VI_
#define _VI_

using arrr =  Eigen::ArrayXd;
typedef Eigen::ArrayXXd arrmat;
typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

#include <boost/math/special_functions/digamma.hpp>
#include "../auxiliary/beta_rng.hpp"

class VI
{

	private:
		const vec &y;
		const mat &X;
		const double tau_a;
		const double tau_b;
		const double omega_a;
		const double omega_b;
		const int max_iter;
		const int N;
		//const int n;
		const int Q;
		double tau_a_n;
		double tau_b_n;
		double omega_a_n;
		double omega_b_n;
		double tau;
		double omega;
		double current_bound;
		double current_iter;
		arrr cluster_count;
		arrr iter_cluster_assignment;
		arrr u;
		vec beta;
		mat V;
	public:
		VI(const vec &y,
		   const mat &X,
		   const double &tau_a,
		   const double &tau_b,
		   const double &omega_a,
		   const double &omega_b,
		   const int &max_iter,
		   std::mt19937 &rng) :
			y(y), X(X), 
			tau_a(tau_a),
			tau_b(tau_b),
			omega_a(omega_a),
			omega_b(omega_b),
			max_iter(max_iter),
			N(y.rows()),
			Q(X.cols())
	{
		initialize_pars(rng);	
	}

		void initialize_pars(std::mt19937 &rng);

		void descend_gradient();

		double calculate_bound();

		bool has_not_converged(arrr &bounds);

		void draw_pars(std::mt19937 &rng,
					   const int &num_samples,
					   arrmat &beta_samples,
					   arrr &sigma_samples,
					   arrr &tau_samples,
					   arrmat &yhat_samples);

		void stick_break(std::mt19937 &rng);

};

#include "VI.inl"


#endif 
