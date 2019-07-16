// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#ifndef heom_eigen_decomposition_hpp
#define heom_eigen_decomposition_hpp

#include "heom/common.hpp"

#include "Eigen/Eigenvalues"

namespace heom {

class eigen_decomposition {
public:
	using eigen_matrix_t = Eigen::MatrixXcd; // NOTE: Eigen defaults to column major format
	//using eigen_map_t = Eigen::Map<eigen_matrix_t>;
	using eigen_map_t = Eigen::Map<Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>; // NOTE: we use row major

	eigen_decomposition(const complex_matrix_t& input);
	complex_vector_t eigenvalues() const; // column vector of eigen values
	complex_matrix_t eigendiag() const; // diagonal matrix of eigen values
	complex_matrix_t eigenvectors() const; // matrix with eigenvectors as columns

	// conversion helpers
	static eigen_matrix_t heom_matrix_to_eigen(complex_matrix_t);
	static complex_matrix_t eigen_matrix_to_heom(eigen_matrix_t);
	static complex_vector_t eigen_vector_to_heom(eigen_matrix_t);

private:
	using eigen_solver_t = Eigen::ComplexEigenSolver<eigen_matrix_t>;

	eigen_solver_t solver_;
};

} // namespace heom

#endif // heom_eigen_decomposition_hpp