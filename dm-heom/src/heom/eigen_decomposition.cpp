// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/eigen_decomposition.hpp"

namespace heom {

eigen_decomposition::eigen_decomposition(const complex_matrix_t& input)
  : solver_(heom_matrix_to_eigen(input), true)
{
}

// column vector of eigen values
complex_vector_t eigen_decomposition::eigenvalues() const
{
	return eigen_vector_to_heom(solver_.eigenvalues());
}

// diagonal matrix of eigenvalues
complex_matrix_t eigen_decomposition::eigendiag() const
{
	return eigen_matrix_to_heom(solver_.eigenvalues().asDiagonal());
}

// matrix with eigenvectors as columns
complex_matrix_t eigen_decomposition::eigenvectors() const
{
	return eigen_matrix_to_heom(solver_.eigenvectors());
}

eigen_decomposition::eigen_matrix_t eigen_decomposition::heom_matrix_to_eigen(complex_matrix_t m)
{
	return eigen_map_t { m.data(), static_cast<Eigen::Index>(m.rows()), static_cast<Eigen::Index>(m.cols()) };
}

complex_matrix_t eigen_decomposition::eigen_matrix_to_heom(eigen_matrix_t m)
{
	complex_matrix_t result { static_cast<size_t>(m.rows()), static_cast<size_t>(m.cols()) };
	for (size_t i = 0; i < result.rows(); ++i) {
		for (size_t j = 0; j < result.cols(); ++j) {
			result.at(i, j) = m(i, j);
		}
	}

	return result;
}

complex_vector_t eigen_decomposition::eigen_vector_to_heom(eigen_matrix_t m)
{
	complex_vector_t result { static_cast<size_t>(m.rows()) };
	for (size_t i = 0; i < result.size(); ++i) {
			result.at(i) = m(i, 0);
	}

	return result;
}


} // namespace heom