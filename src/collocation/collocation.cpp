// 
//! \file collocation.cpp
//  src
//  
//  Created by Justin Brown on 2013-04-02.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

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
		LOG4CXX_TRACE (config::logger, "Indexing...");
		
		LOG4CXX_TRACE (config::logger, "Value is " << data [deriv * rows * cols + row * cols + col])
		
		return data [deriv * rows * cols + row * cols + col];
	}

	cheb_grid::cheb_grid (int i_M, int i_N) : collocation_grid (3, i_M, i_N) {
		int d, m, k;
		pioN = std::acos (-1.0) / i_N;
		
		LOG4CXX_TRACE (config::logger, "Instantiating...");
		
		for (d = 0; d < 3; ++d) {
			for (m = 0; m < i_M; ++m) {
				for (k = 0; k < i_N; ++k) {
					LOG4CXX_TRACE (config::logger, "" << d << " " << m << " " << k);
					LOG4CXX_TRACE (config::logger, "" << index (d, m, k));
					index (d, m, k) = recursion (d, m, k);
					exists (d, m, k) = true;
				}
			}
		}
		
		LOG4CXX_TRACE (config::logger, "Instantiated...")
	}
	
	double cheb_grid::recursion (int d, int m, int k) {
		LOG4CXX_TRACE (config::logger, "Calculating d^" << d << " T_" << m << " (x_" << k << ")...");

		if (exists (d, m, k)) {
			LOG4CXX_DEBUG (config::logger, "Already exists, using table value...");
			return index (d, m, k);
		} else if ((d == 2 && (m == 0 || m == 1) || (d == 1 && m == 0))) {
			LOG4CXX_DEBUG (config::logger, "It's 0.0...");
			return 0.0;
		} else if ((d == 1 && m == 1) || (d == 0 && m == 0)) {
			LOG4CXX_DEBUG (config::logger, "First one...");
			return 1.0;
		} else if (d == 0 && m == 1) {
			LOG4CXX_DEBUG (config::logger, "It's x...");
			return pioN * k;
		} else if (d == 0) {
			LOG4CXX_DEBUG (config::logger, "Recursing cheb...");
			return 2.0 * pioN * k * recursion (0, m - 1, k) + recursion (0, m - 2, k);
		} else if (d == 1) {
			LOG4CXX_DEBUG (config::logger, "Recursing dcheb...");
			return 2.0 * recursion (0, m - 1, k) + 2.0 * pioN * k * recursion (1, m - 1, k) - recursion (1, m - 2, k);
		} else if (d == 2) {
			LOG4CXX_DEBUG (config::logger, "Recursing d2cheb...");
			return 4.0 * recursion (1, m - 1, k) + 2.0 * pioN * k * recursion (2, m - 1, k) - recursion (2, m - 2, k);
		} else {
			LOG4CXX_DEBUG (config::logger, "Well, damn...");
			throw 0;
		}
		
		LOG4CXX_TRACE (config::logger, "Calculated.");
	}
} /* collocation */