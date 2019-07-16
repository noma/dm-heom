// This file is part of DM-HEOM (https://github.com/noma/dm-heom)
//
// Copyright (c) 2015-2019 Matthias Noack, Zuse Institute Berlin
// Copyright (c) 2015 Tobias Kramer, Zuse Institute Berlin
//
// Licensed under the 3-clause BSD License, see accompanying LICENSE,
// CONTRIBUTORS.md, and README.md for further information.

#include "heom/integer_partition.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace heom {
namespace ip {

// recursive definition
// http://fr.wikipedia.org/wiki/Partition_d%27un_entier
int_t size_partiton_restrict(int_t n, int_t k)
{
	if ((n == 0) && (k == 0)) return 1;
	else if ((n <= 0) || (k <= 0)) return 0;
	else return size_partiton_restrict(n - k, k) + size_partiton_restrict(n - 1, k - 1);
}

// summing up all the partition sizes to obtain the total number of
// restricted partitions
int_t size_my_integer_partitions(int_t n, int_t s)
{
	// compute the total number of restricted partitions
	int_t res = 0;
	for (int_t p = 1; p <= s; ++p) {
		res += size_partiton_restrict(n, p);
	}
	return res;
}

// generate all integer partitions of n with maximum number of elements s
// returns the number of generated partitions
int_t my_integer_partitions(int_t n, int_t s, int_t* res)
{
	int_t i, r, h, m, p, j, resi; //, max
	std::unique_ptr<int_t[]> x{ new int_t[n + 1] };

	const int_t SIZE = s;

	// compute the expected number of partitions
	int_t resimax = size_my_integer_partitions(n, s);

	resi = 0;

	for (i = 1; i <= n; ++i)
		x[i] = 1;

	if (n <= s) {
		for (p = 0; p < n; ++p) res[resi * SIZE + p] = x[p + 1];
		for (; p < s; ++p) res[resi * SIZE + p] = 0;
		++resi;
	}

	if (n > 1) {
		x[0] = -1;
		x[1] = 2;
		h = 1;
		m = n - 1;
		if (m <= s) {
			for (p = 0; p < m; ++p) res[resi * SIZE + p] = x[p + 1];
			for (; p < s; ++p) res[resi * SIZE + p] = 0;
			++resi;

		}
		while (x[1] != n) {
			if (m - h > 1) {
				h = h + 1;
				x[h] = 2;
				m = m - 1;
			} else {
				j = m - 2;
				while (x[j] == x[m - 1]) {
					x[j] = 1;
					j = j - 1;
				}
				h = j + 1;
				x[h] = x[m - 1] + 1;
				r = x[m] + x[m - 1] * (m - h - 1);
				x[m] = 1;
				if (m - h > 1) {
					x[m - 1] = 1;
				}
				m = h + r - 1;
			}

			if (m <= s) {
				for (p = 0; p < m; ++p) res[resi * SIZE + p] = x[p + 1];
				for (; p < s; ++p) res[resi * SIZE + p] = 0;
				++resi;
			}
		}
	}
	if (resi != resimax) {
		throw std::range_error("ip::my_integer_partitions(): error: wrong number of partitions.");
	}

	return resi;
}

} // namespace ip
} // namesapce heom
