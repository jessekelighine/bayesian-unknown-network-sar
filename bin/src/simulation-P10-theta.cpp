// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Hamiltonian.h>
#include <Exchange.h>
#include <ERGM_stat.h>
#include <vector>
#include <cmath>
#include <iostream>

namespace Hamiltonian
{
	typedef std::pair<int,arma::vec> ParamsPair;

	arma::umat ToZeroDiagonalBinaryMatrix (
			const arma::ivec& discrete_positions,
			const int size
			) {
		arma::umat output (size, size, arma::fill::zeros);
		for ( int j { 0 }; j < std::pow(size,2) ; j++ ) {
			if ( j / size == j % size ) continue;
			if (  discrete_positions(j) < 0 ) continue;
			output(j) = 1;
		}
		return output;
	}

	template<>
	double PotentialEnergyDifference<ParamsPair> (
			const arma::ivec& discrete_positions,
			const int coord,
			const ParamsPair& params
			) {
		const int coord_row = coord % params.first;
		const int coord_col = coord / params.first;
		if ( coord_row == coord_col ) return 0;
		const arma::vec stats {
			static_cast<double>(ERGM_stat::diff_edges(coord, discrete_positions)),
			static_cast<double>(ERGM_stat::diff_mutual(coord, discrete_positions, params.first))
		};
		return arma::dot(stats, params.second);
	}

}

namespace Exchange
{
	template<>
	double LogLikelihood<arma::umat> (
			const arma::umat& data,
			const arma::vec& parameters
			) {
		arma::vec stats {
			static_cast<double>(ERGM_stat::edges(data)),
			static_cast<double>(ERGM_stat::mutual(data))
		};
		return arma::dot(stats, parameters);
	}

	template<>
	arma::umat AuxSample<arma::umat> (
			const arma::vec& parameters,
			const int size
			) {
		const int dimension = std::pow(size, 2);
		const int chain_length = 5;
		const int cycle_per_iteration = 5;
		const Hamiltonian::ParamsPair& params = { size, parameters };
		/* Hamiltonian Run */
		arma::vec positions = Hamiltonian::MomentumSampler(dimension);
		for ( int iteration { 0 }; iteration < chain_length; iteration++ ) {
			Hamiltonian::PositionDynamics(positions, dimension, params, cycle_per_iteration);
		}
		return Hamiltonian::ToZeroDiagonalBinaryMatrix(
				Hamiltonian::ToDiscrete(positions, dimension),
				params.first);
	}
}

// [[Rcpp::export]]
arma::vec SamplePosteriorTheta (
		const arma::umat& data,
		const int parameters_dimension,
		const arma::vec& prior_mean,
		const double prior_variance = 1,
		const int chain_length = 800
		) {
	arma::vec state (parameters_dimension, arma::fill::zeros);
	for ( int iteration { 1 }; iteration < chain_length ; iteration++ ) { 
		Exchange::RunOnce(Exchange::Propose(state, 0.25), state, data, prior_mean, prior_variance);
	}
	return state;
}
