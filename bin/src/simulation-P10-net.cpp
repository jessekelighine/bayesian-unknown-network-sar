// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <LogLikelihood.h>
#include <Hamiltonian.h>
#include <ERGM_stat.h>

namespace LogLikelihood
{
	double KernelAll (
			const Rcpp::List data,
			const double lambda,
			const arma::mat& W,
			const arma::vec& beta1,
			const arma::vec& beta2,
			const double sigma2
			) {
		const Rcpp::List& y    = data["Y"];
		const Rcpp::List& X1   = data["X1"];
		const Rcpp::List& X2   = data["X2"];
		const int panel_length = data["panel_length"];
		double output = 0;
		for ( int t { 0 }; t < panel_length; t++ ) {
			const arma::vec& y_t  = y[t];
			const arma::mat& X1_t = X1[t];
			const arma::mat& X2_t = X2[t];
			output += LogLikelihood::Kernel(y_t, X1_t, X2_t, lambda, W, beta1, beta2, sigma2);
		}
		return output;
	}
}

namespace Hamiltonian
{
	template<typename T>
	T ToBinaryMatrix (
			const arma::ivec& discrete_positions,
			const int size
			) {
		T output (size, size, arma::fill::zeros);
		for ( int j { 0 }; j < std::pow(size,2); j++ ) {
			if ( j / size == j % size ) continue;
			if ( discrete_positions(j) < 0 ) continue;
			output(j) = 1;
		}
		return output;
	}

	arma::mat ToNormalizedMatrix (
			const arma::ivec& discrete_positions,
			const int size,
			const int coord = -1
			) {
		arma::mat output = ToBinaryMatrix<arma::mat>(discrete_positions, size);
		if ( coord >= 0 ) output(coord) = ( output(coord) > 0 ) ? 0 : 1;
		return arma::normalise(output, 1, 1);
	}

	template<>
	double PotentialEnergyDifference<Rcpp::List> (
			const arma::ivec& discrete_positions,
			const int coord,
			const Rcpp::List& params
			) {
		const Rcpp::List& data  = params["data"];
		const int size = data["size"];
		const int coord_row = coord % size;
		const int coord_col = coord / size;
		if ( coord_col == coord_row ) return 0;
		/* Extract Params */
		const arma::vec& theta = params["theta"];
		const arma::vec& beta1 = params["beta1"];
		const arma::vec& beta2 = params["beta2"];
		const double lambda = params["lambda"];
		const double sigma2 = params["sigma2"];
		/* SAR Matrix */
		const arma::mat W_s = ToNormalizedMatrix(discrete_positions, size);
		const arma::mat W_p = ToNormalizedMatrix(discrete_positions, size, coord);
		/* Log Determinant */
		double logdet_W_s, logdet_W_p, sig;
		arma::log_det( logdet_W_s, sig, LogLikelihood::CalcM(lambda,W_s) );
		arma::log_det( logdet_W_p, sig, LogLikelihood::CalcM(lambda,W_p) );
		/* Sufficient Statistics */
		const arma::vec sufficient_stats {
			static_cast<double>(ERGM_stat::diff_edges(coord, discrete_positions)),
			static_cast<double>(ERGM_stat::diff_mutual(coord, discrete_positions, size))
		};
		/* Calculate Output */
		double output = 0;
		output += logdet_W_p - logdet_W_s;
		output += LogLikelihood::KernelAll(data, lambda, W_p, beta1, beta2, sigma2);
		output -= LogLikelihood::KernelAll(data, lambda, W_s, beta1, beta2, sigma2);
		output += arma::dot(sufficient_stats, theta);
		return output;
	}
}

// [[Rcpp::export]]
arma::umat SamplePosteriorNet (
		const int dimension,
		const Rcpp::List& params,
		const int chain_length = 1000,
		const int cycle_per_iteration = 10
		) {
	arma::vec positions = Hamiltonian::MomentumSampler(dimension);
	for ( int j { 0 }; j < chain_length; j++ ) {
		Hamiltonian::PositionDynamics(positions, dimension, params, cycle_per_iteration);
	}
	return Hamiltonian::ToBinaryMatrix<arma::umat>(
			Hamiltonian::ToDiscrete(positions, dimension),
			static_cast<int>(std::sqrt(dimension))
			);
}
