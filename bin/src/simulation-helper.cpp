// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat RowNormalize (
		arma::mat M,
		const bool zero_diag = true
		) {
	if ( zero_diag ) M.diag().zeros();
	return arma::normalise(M, 1, 1);
}
