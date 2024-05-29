#ifndef EXCHANGE_H
#define EXCHANGE_H

#include <RcppArmadillo.h>
#include <iostream>

namespace Exchange
{
	/* TO BE DEFINED BY USER 
	 *
	 * This function returns the log-likleihood of the 'data' given the
	 * 'parameters'.
	 *
	 * */
	template<typename T>
	double LogLikelihood (
			const T& data,
			const arma::vec& parameters
			);

	/* TO BE DEFINED BY USER
	 *
	 * This function returns a sample from of "data" given the
	 * 'parameters'.
	 *
	 * */
	template<typename T>
	T AuxSample (
			const arma::vec& parameters,
			const int size
			);

	double PriorLogDensity (
			const arma::vec& parameters,
			const arma::vec& prior_mean,
			const double prior_variance = 30
			) {
		double output = 0;
		for ( int i { 0 }; i < parameters.n_rows; i++ ) {
			// output += -std::log(2 * M_PI)/2
			output += -std::log(prior_variance) / 2;
			output += -std::pow(parameters(i) - prior_mean(i), 2) / ( 2 * prior_variance );
		}
		return output;
	}

	arma::vec Propose (
			const arma::vec& parameters,
			const double sd = 0.025
			) {
		arma::vec output (parameters.n_rows, arma::fill::randn);
		output = output * sd + parameters;
		return output;
	}

	template<typename T>
	void RunOnce (
			const arma::vec& proposal,
			arma::vec& state,
			const T& data,
			const arma::vec& prior_mean,
			const double prior_variance = 30
			) {
		const T aux_data = AuxSample<T>(proposal, data.n_rows);
		double metropolis_ratio = 0;
		metropolis_ratio += LogLikelihood(data, proposal);
		metropolis_ratio -= LogLikelihood(data, state);
		metropolis_ratio += PriorLogDensity(proposal, prior_mean, prior_variance);
		metropolis_ratio -= PriorLogDensity(state, prior_mean, prior_variance);
		metropolis_ratio += LogLikelihood(aux_data, state);
		metropolis_ratio -= LogLikelihood(aux_data, proposal);
		if ( std::log(arma::randu()) < metropolis_ratio ) state = proposal;
	}
}

#endif
