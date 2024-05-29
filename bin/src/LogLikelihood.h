#ifndef LOGLIKELIHOOD_H
#define LOGLIKELIHOOD_H

#include <cmath>
#include <vector>
#include <RcppArmadillo.h>

namespace LogLikelihood
{
	arma::mat CalcM (
			const double lambda,
			const arma::mat& W
		     ) {
		return arma::eye(arma::size(W)) - lambda * W;
	}

	double Kernel (
			const arma::vec& y,
			const arma::mat& X1,
			const arma::mat& X2,
			const double lambda,
			const arma::mat& W,
			const arma::vec& beta1,
			const arma::vec& beta2,
			const double sigma2
			) {
		const arma::vec reduced = CalcM(lambda,W) * y - X1 * beta1 - W * X2 * beta2;
		return -arma::dot(reduced, reduced) / ( 2 * sigma2 );
	}

	double Group (
			const arma::vec& y,
			const arma::mat& X1,
			const arma::mat& X2,
			const double lambda,
			const arma::mat& W,
			const arma::vec& beta1,
			const arma::vec& beta2,
			const double sigma2
			) {
		// Same as Kernel
		const int group_size = y.n_rows;
		const arma::mat M = CalcM(lambda, W);
		const arma::vec reduced = M * y - X1 * beta1 - W * X2 * beta2;
		double kernel = -arma::dot(reduced, reduced) / ( 2 * sigma2 );
		// Determinant
		double log_det_M, sig;
		arma::log_det(log_det_M, sig, M);
		// Output
		return log_det_M - group_size * std::log(sigma2) / 2 + kernel;
	}
}

#endif
