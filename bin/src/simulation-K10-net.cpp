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
	typedef std::pair<int,int> PairRowCol;

	PairRowCol Coord2RowCol (
			const int coord,
			const int size
			) {
		const int i = size - 2 - std::floor(std::sqrt(-8 * coord + 4 * size * (size - 1) - 7) / 2.0 - 0.5);
		const int j = coord + i + 1 - size * (size - 1) / 2 + (size - i) * (size - i - 1) / 2;
		return PairRowCol (j,i);
	}

	// template<typename T>
	// T ToBinaryMatrix (
	// 		const arma::ivec& discrete_positions,
	// 		const int size
	// 		) {
	// 	T output (size, size, arma::fill::zeros);
	// 	for ( int j { 0 }; j < std::pow(size,2); j++ ) {
	// 		if ( j / size >= j % size ) continue;
	// 		if ( discrete_positions(coord_iter) > 0 ) output(coord_iter) = 1;
	// 		coord_iter++;
	// 	}
	// 	return output + output.t();
	// }

	template<typename T>
	T ToBinaryMatrix (
			const arma::ivec& discrete_positions,
			const int size,
			const int coord = -1
			) {
		T output (size, size, arma::fill::zeros);
		for ( int j { 0 }; j < (size * size - size) / 2; j++ ) {
			if ( discrete_positions(j) < 0 ) continue;
			PairRowCol coord_rowcol = Coord2RowCol(j, size);
			output(coord_rowcol.first, coord_rowcol.second) = 1;
		}
		return output + output.t();
	}

	arma::umat FlipCoord (
			arma::umat matrix,
			const int coord,
			const int size
			) {
		const PairRowCol coord_rowcol = Coord2RowCol(coord, size);
		const int change = ( matrix(coord_rowcol.first, coord_rowcol.second) > 0 ) ? 0 : 1;
		matrix(coord_rowcol.first, coord_rowcol.second) = change;
		matrix(coord_rowcol.second, coord_rowcol.first) = change;
		return matrix;
	}

	template<>
	double PotentialEnergyDifference<Rcpp::List> (
			const arma::ivec& discrete_positions,
			const int coord,
			const Rcpp::List& params
			) {
		/* Extract Params */
		const Rcpp::List& data = params["data"];
		const arma::vec& theta = params["theta"];
		const arma::vec& beta1 = params["beta1"];
		const arma::vec& beta2 = params["beta2"];
		const double lambda = params["lambda"];
		const double sigma2 = params["sigma2"];
		const int size = data["size"];
		/* SAR Matrix */
		const arma::umat net_s = ToBinaryMatrix<arma::umat>(discrete_positions, size);
		const arma::umat net_p = FlipCoord(net_s, coord, size);
		/* Sufficient Statistics */
		const arma::vec sufficient_stats {
			static_cast<double>(ERGM_stat::edges(net_p)    - ERGM_stat::edges(net_s)),
			static_cast<double>(ERGM_stat::kstar(net_p,2)  - ERGM_stat::kstar(net_s,2)),
			static_cast<double>(ERGM_stat::kstar(net_p,3)  - ERGM_stat::kstar(net_s,3)),
			static_cast<double>(ERGM_stat::triangle(net_p) - ERGM_stat::triangle(net_s))
		};
		/* SAR Matrix */
		const arma::mat W_s = arma::normalise(arma::conv_to<arma::mat>::from(net_s), 1, 1);
		const arma::mat W_p = arma::normalise(arma::conv_to<arma::mat>::from(net_p), 1, 1);
		/* Log Determinant */
		double logdet_W_s, logdet_W_p, sig;
		arma::log_det( logdet_W_s, sig, LogLikelihood::CalcM(lambda,W_s) );
		arma::log_det( logdet_W_p, sig, LogLikelihood::CalcM(lambda,W_p) );
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
			Rcpp::as<int>(params["size"])
			);
}
