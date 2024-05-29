#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <cmath>
#include <vector>
#include <RcppArmadillo.h>

namespace Hamiltonian
{
	typedef std::pair<int,double> CoordHitTimePair;

	arma::vec MomentumSampler (
			const int dimension
			) {
		arma::vec output (dimension);
		output.randn(dimension);
		return output;
	}

	double MomentumDynamics (
			const int coord,
			const double run_time,
			const arma::vec& positions,
			const arma::vec& momentums
			) {
		return momentums(coord) * std::cos(run_time) - positions(coord) * std::sin(run_time);
	}

	double CalcHitTime (
			const int coord,
			const arma::vec& positions,
			const arma::vec& momentums
			) {
		const double phi = std::atan2(positions(coord), momentums(coord));
		if ( phi > 0 ) return M_PI - phi;
		return -phi;
	}

	arma::ivec ToDiscrete (
			const arma::vec& positions,
			const int dimension
			) {
		arma::ivec output (dimension);
		for ( int coord { 0 }; coord < dimension ; coord++ ) {
			output(coord) = ( positions(coord) > 0 ) ? 1 : -1;
		}
		return output;
	}

	/* TO BE DEFINED BY USER
	 *
	 * This function should return the difference in "log-likelihood" when
	 * the 'coord'-th entry of 'discrete_positions' is flipped. That is,
	 * let $L(p)$ denote the log-likleihood of 'discrete_positions' $p$ and
	 * let $p'$ denote the 'discrete_positions' with the 'coord'-th entry
	 * flipped, then this function should return $L(p') - L(p)$.
	 *
	 * */
	template<typename T>
		double PotentialEnergyDifference (
				const arma::ivec& discrete_positions,
				const int coord,
				const T& params
				);

	template<typename T>
		double CalcMomentumPlusSquared (
				const arma::vec& momentums,
				const arma::ivec& discrete_positions,
				const int coord,
				const T& params
				) {
			double output = 0;
			output += 2 * PotentialEnergyDifference<T>(discrete_positions, coord, params);
			output += pow(momentums(coord), 2);
			return output;
		}

	template<typename T>
		void PositionDynamics (
				arma::vec& positions,
				const int dimension,
				const T& params,
				const int total_cycles
				) {
			arma::vec momentums = MomentumSampler(dimension);
			/* Calculate hit times and order coordinates by hit time */
			std::vector<CoordHitTimePair> coord_hittime_pairs (dimension);
			for ( int coord { 0 }; coord < dimension ; coord++ ) {
				coord_hittime_pairs[coord] = std::make_pair(coord, CalcHitTime(coord, positions, momentums));
			}
			std::sort(
					coord_hittime_pairs.begin(),
					coord_hittime_pairs.end(),
					[] ( const CoordHitTimePair& pair1, const CoordHitTimePair& pair2 ) { return pair1.second < pair2.second; }
				 );
			/* Binary Position */
			arma::ivec discrete_positions = ToDiscrete(positions, dimension);
			/* First Cycle */
			for ( const auto& pair : coord_hittime_pairs ) {
				momentums(pair.first) = MomentumDynamics( pair.first, pair.second, positions, momentums);
				const double momentum_plus_squared = CalcMomentumPlusSquared( momentums, discrete_positions, pair.first, params);
				if ( momentum_plus_squared > 0 ) {
					discrete_positions(pair.first) *= -1;
					momentums(pair.first) = discrete_positions(pair.first) * std::sqrt(momentum_plus_squared);
				} else {
					momentums(pair.first) *= -1;
				}
			}
			/* 2 to total_cycles-1 cycles */
			for ( int cycle { 1 }; cycle < total_cycles ; cycle++ ) {
				for ( const auto& pair : coord_hittime_pairs ) {
					const double momentum_plus_squared = CalcMomentumPlusSquared(momentums, discrete_positions, pair.first, params);
					if ( momentum_plus_squared > 0 ) {
						discrete_positions(pair.first) *= -1;
						momentums(pair.first) = discrete_positions(pair.first) * std::sqrt(momentum_plus_squared);
					}
				}
			}
			/* Last Cycle */
			for ( const auto& pair : coord_hittime_pairs ) {
				double final_run_time;
				if ( pair.second < M_PI/2 ) {
					const double momentum_plus_squared = CalcMomentumPlusSquared(momentums, discrete_positions, pair.first, params);
					if ( momentum_plus_squared > 0 ) {
						discrete_positions(pair.first) *= -1;
						momentums(pair.first) = discrete_positions(pair.first) * std::sqrt(momentum_plus_squared);
					}
					final_run_time = M_PI * 1/2 - pair.second;
				} else {
					final_run_time = M_PI * 3/2 - pair.second;
				}
				positions(pair.first) = momentums(pair.first) * std::sin(final_run_time);
			}
		}
}

#endif
