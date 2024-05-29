#ifndef ERGM_STAT_H
#define ERGM_STAT_H

#include <RcppArmadillo.h>

namespace ERGM_stat
{
	int choose ( const int n, const int k ) { return ( k == 0 ) ? 1 : choose(n-1, k-1) * n / k; }

	/* Directed Networks */

	int edges  ( const arma::umat& matrix ) { return arma::accu(matrix); }
	int mutual ( const arma::umat& matrix ) { return arma::accu((matrix == matrix.t()) && (matrix || matrix.t())) / 2; }
	int ttriad ( const arma::umat& matrix ) { return arma::trace(matrix * matrix * matrix.t()); }
	int ctriad ( const arma::umat& matrix ) { return arma::trace(matrix * matrix * matrix) / 3; }

	/* Undirected Networks */

	int triangle ( const arma::umat& matrix ) { return arma::trace(matrix * matrix * matrix) / 6; }
	int kstar ( const arma::umat& matrix, const int k ) {
		arma::uvec diag = arma::diagvec(matrix * matrix);
		return arma::accu(diag.transform( [k] ( const int i ) { return choose(i,k); } ));
	}

	/* Undirected Networks */

	int diff_kstar (
			const int i,
			const int j,
			const arma::umat& matrix,
			const int k_kstar
			) {
		const int diag_i  = arma::accu(matrix.col(i));
		const int diag_j  = arma::accu(matrix.col(j));
		const int diag_i_ = diag_i + 2 * ( ( matrix(i,j) == 0 ) ? 0 : -1 ) + 1;
		const int diag_j_ = diag_j + 2 * ( ( matrix(i,j) == 0 ) ? 0 : -1 ) + 1;
		int output = 0;
		output += ( diag_i_ >= k_kstar ) ? choose(diag_i_, k_kstar) : 0;
		output += ( diag_j_ >= k_kstar ) ? choose(diag_j_, k_kstar) : 0;
		output -= ( diag_i  >= k_kstar ) ? choose(diag_i,  k_kstar) : 0;
		output -= ( diag_j  >= k_kstar ) ? choose(diag_j,  k_kstar) : 0;
		return output;
	}

	int diff_triangle (
			const int i,
			const int j,
			const arma::umat& matrix
			) {
		const int output = arma::dot(matrix.col(i), matrix.col(j));
		return ( matrix(i,j) == 0 ) ? output : -output;
	}

	int diff_edges (
			const int i,
			const int j,
			const arma::umat& matrix
			) {
		return ( matrix(i,j) == 0 ) ? 1 : -1;
	}

	/* Directed Networks */

	int diff_edges (
			const int coord,
			const arma::ivec& discrete_positions
			) {
		return -discrete_positions(coord);
	}

	int diff_mutual (
			const int coord,
			const arma::ivec& discrete_positions,
			const int size
			) {
		const int t_coord = ( coord / size ) + ( coord % size ) * size;
		return -discrete_positions(coord) * static_cast<int>(discrete_positions(t_coord)==1);
	}

	int diff_ttriad (
			const int coord,
			const arma::ivec& discrete_positions,
			const int size
			) {
		int output = 0;
		const int i = coord % size;
		const int j = coord / size;
		if ( i == j ) return 0;
		for ( int k { 0 }; k < size; k++ ) {
			const int a_ik = ( discrete_positions(i + k * size) > 0  && i != k ) ? 1 : 0;
			const int a_kj = ( discrete_positions(k + j * size) > 0  && k != j ) ? 1 : 0;
			const int a_jk = ( discrete_positions(j + k * size) > 0  && j != k ) ? 1 : 0;
			const int a_ki = ( discrete_positions(k + i * size) > 0  && k != i ) ? 1 : 0;
			output += a_ik * a_kj + a_ki * a_kj + a_ik * a_jk;
		}
		return ( discrete_positions(coord) > 0 ) ? -output : output;
	}

	int diff_ctriad (
			const int coord,
			const arma::ivec& discrete_positions,
			const int size
			) {
		int output = 0;
		const int i = coord % size;
		const int j = coord / size;
		if ( i == j ) return 0;
		for ( int k { 0 }; k < size; k++ ) {
			const int a_jk = ( discrete_positions(j + k * size) > 0 && j != k ) ? 1 : 0;
			const int a_ki = ( discrete_positions(k + i * size) > 0 && k != i ) ? 1 : 0;
			output += a_jk * a_ki;
		}
		return ( discrete_positions(coord) > 0 ) ? -output : output;
	}
}

#endif
