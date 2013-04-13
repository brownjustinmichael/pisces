/*!***********************************************************************
 * \file collocation.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cassert>
#include <exception>
#include "../config.hpp"
#include "collocation.hpp"

namespace collocation
{
	collocation_grid::collocation_grid (int i_derivs, int i_rows, int i_cols) {
		int i;
		rows = i_rows;
		cols = i_cols;
		derivs = i_derivs;
		
		TRACE ("Instantiating...")
		
		data.resize (derivs);
		
		for (i = 0; i < i_derivs; ++i) {
			data [i].resize (i_rows * i_cols);
		}
		
		TRACE ("Instantiated...")
	}

	chebyshev_grid::chebyshev_grid (int i_M, int i_N, double i_scale) : collocation_grid (3, i_M, i_N) {
		int d, m, k;
		scale = i_scale;
		pioN = std::acos (-1.0) / i_N;
		exists_array.resize (i_M * i_N * 3, false);
		
		TRACE ("Instantiating...");
		
		for (d = 0; d < 3; ++d) {
			for (m = 0; m < i_M; ++m) {
				for (k = 0; k < i_N; ++k) {
					index (d, m, k) = recursion (d, m, k);
					exists (d, m, k) = true;
				}
			}
		}
		
		for (d = 0; d < 3; ++d) {
			for (k = 0; k < i_N; ++k) {
				index (d, 0, k) /= 2.0;
				index (d, i_M - 1, k) /= 2.0;
			}
		}
		
		TRACE ("Instantiated...")
	}
	
	double chebyshev_grid::recursion (int d, int m, int k) {		
		assert (d >= 0);
		assert (d < 3);
		assert (m >= 0);
		assert (k >= 0);

		// Use the recursion relations to calculate the correct value of the polynomial at the collocation point
		if (exists (d, m, k)) {
			// This value has already been calculated
			return index (d, m, k);
		} else if ((d == 2 && (m == 0 || m == 1)) || (d == 1 && m == 0)) {
			// The first two polynomials of the second derivative and the first polynomial of the first derivative are 0.0
			return 0.0;
		} else if ((d == 1 && m == 1) || (d == 0 && m == 0)) {
			// The second polynomial of the first derivative and the first Chebyshev polynomial are 1.0
			return scale;
		} else if (d == 0 && m == 1) {
			// The second Chebyshev polynomial is cos (theta)
			return scale * std::cos (pioN * k);
		} else if (d == 0) {
			// Use recursion to find the Chebyshev polynomial
			return 2.0 * std::cos (pioN * k) * recursion (0, m - 1, k) - recursion (0, m - 2, k);
		} else if (d == 1) {
			// Use recursion to find the first derivative of the Chebyshev polynomial
			return 2.0 * recursion (0, m - 1, k) + 2.0 * std::cos (pioN * k) * recursion (1, m - 1, k) - recursion (1, m - 2, k);
		} else if (d == 2) {
			// Use recursion to find the second derivative of the Chebyshev polynomial
			return 4.0 * recursion (1, m - 1, k) + 2.0 * std::cos (pioN * k) * recursion (2, m - 1, k) - recursion (2, m - 2, k);
		} else {
			// There's been a bad index, kill the program
			FATAL ("Bad index")
			throw 0;
		}
	}
} /* collocation */