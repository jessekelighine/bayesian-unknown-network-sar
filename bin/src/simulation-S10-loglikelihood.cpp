// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <LogLikelihood.h>

// [[Rcpp::export]]
double LogLikelihoodAll (
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
		output += LogLikelihood::Group(y_t, X1_t, X2_t, lambda, W, beta1, beta2, sigma2);
	}
	return output;
}

// [[Rcpp::export]]
double LogLikelihoodKernelAll (
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
