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
	typedef std::pair<int,int> PairRowCol;

	PairRowCol Coord2RowCol (
			const int coord,
			const int size
			) {
		const int i = size - 2 - std::floor(std::sqrt(-8 * coord + 4 * size * (size - 1) - 7) / 2.0 - 0.5);
		const int j = coord + i + 1 - size * (size - 1) / 2 + (size - i) * (size - i - 1) / 2;
		return PairRowCol (j,i);
	}

	arma::umat ToBinaryMatrix (
			const arma::ivec& discrete_positions,
			const int size
			) {
		arma::umat output (size, size, arma::fill::zeros);
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
	double PotentialEnergyDifference<ParamsPair> (
			const arma::ivec& discrete_positions,
			const int coord,
			const ParamsPair& params
			) {
		const arma::umat net_s = ToBinaryMatrix(discrete_positions, params.first);
		const arma::umat net_p = FlipCoord(net_s, coord, params.first);
		const arma::vec stats {
			static_cast<double>(ERGM_stat::edges(net_p)    - ERGM_stat::edges(net_s)),
			static_cast<double>(ERGM_stat::kstar(net_p,2)  - ERGM_stat::kstar(net_s,2)),
			static_cast<double>(ERGM_stat::kstar(net_p,3)  - ERGM_stat::kstar(net_s,3)),
			static_cast<double>(ERGM_stat::triangle(net_p) - ERGM_stat::triangle(net_s))
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
			static_cast<double>(ERGM_stat::kstar(data,2)),
			static_cast<double>(ERGM_stat::kstar(data,3)),
			static_cast<double>(ERGM_stat::triangle(data))
		};
		return arma::dot(stats, parameters);
	}

	template<>
	arma::umat AuxSample<arma::umat> (
			const arma::vec& parameters,
			const int size
			) {
		const int dimension = ( size * size - size ) / 2;
		const int chain_length = 8;
		const int cycle_per_iteration = 6;
		const Hamiltonian::ParamsPair& params = { size, parameters };
		/* Hamiltonian Run */
		arma::vec positions = Hamiltonian::MomentumSampler(dimension);
		for ( int iteration { 0 }; iteration < chain_length; iteration++ ) {
			Hamiltonian::PositionDynamics(positions, dimension, params, cycle_per_iteration);
		}
		return Hamiltonian::ToBinaryMatrix(
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
		const int chain_length = 1000
		) {
	// arma::vec state (parameters_dimension, arma::fill::zeros); # Works for K1, K2, and K4
	arma::vec state = prior_mean;
	for ( int iteration { 1 }; iteration < chain_length ; iteration++ ) { 
		Exchange::RunOnce(Exchange::Propose(state, 0.25), state, data, prior_mean, prior_variance);
	}
	return state;
}
