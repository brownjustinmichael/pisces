// 
//! \file collocation.cpp
//  src
//  
//  Created by Justin Brown on 2013-04-02.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#include <cassert>
#include <exception>
#include "../config.hpp"
#include "collocation.hpp"

namespace collocation
{
	collocation_grid::collocation_grid (int i_derivs, int i_rows, int i_cols) {
		rows = i_rows;
		cols = i_cols;
		derivs = i_derivs;
		
		LOG4CXX_TRACE (config::logger, "Instantiating...")
		
		data.resize (i_rows * i_cols * i_derivs);
		exists_array.resize (i_rows * i_cols * i_derivs, false);
		
		LOG4CXX_TRACE (config::logger, "Instantiated...")
	}
	
	double &collocation_grid::index (int deriv, int row, int col) {
		return data [deriv * rows * cols + row * cols + col];
	}

	cheb_grid::cheb_grid (int i_M, int i_N) : collocation_grid (3, i_M, i_N) {
		int d, m, k;
		pioN = std::acos (-1.0) / i_N;
		
		LOG4CXX_TRACE (config::logger, "Instantiating...");
		
		for (d = 0; d < 3; ++d) {
			for (m = 0; m < i_M; ++m) {
				for (k = 0; k < i_N; ++k) {
					index (d, m, k) = recursion (d, m, k);
					exists (d, m, k) = true;
				}
			}
		}
		
		LOG4CXX_TRACE (config::logger, "Instantiated...")
	}
	
	double cheb_grid::recursion (int d, int m, int k) {
		LOG4CXX_TRACE (config::logger, "Calculating d^" << d << " T_" << m << " (x_" << k << ")...");
		
		assert (d >= 0);
		assert (d < 3);
		assert (m >= 0);
		assert (k >= 0);

		// Use the recursion relations to calculate the correct value of the polynomial at the collocation point
		if (exists (d, m, k)) {
			// This value has already been calculated
			return index (d, m, k);
		} else if ((d == 2 && (m == 0 || m == 1) || (d == 1 && m == 0))) {
			// The first two polynomials of the second derivative and the first polynomial of the first derivative are 0.0
			return 0.0;
		} else if ((d == 1 && m == 1) || (d == 0 && m == 0)) {
			// The second polynomial of the first derivative and the first Chebyshev polynomial are 1.0
			return 1.0;
		} else if (d == 0 && m == 1) {
			// The second Chebyshev polynomial is cos (theta)
			return std::cos (pioN * k);
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
			LOG4CXX_FATAL (config::logger, "Bad index")
			throw 0;
		}
		
		LOG4CXX_TRACE (config::logger, "Calculated.");
	}
} /* collocation */