// 
//! \file collocation.hpp
//  collocation
//  
//  Created by Justin Brown on 2013-04-02.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef COLLOCATION_HPP_HLJUSAEZ
#define COLLOCATION_HPP_HLJUSAEZ

#include <vector>
#include <cmath>

namespace collocation
{
	class collocation_grid
	{
	private:
		int rows;
		int cols;
		int derivs;
		std::vector<double> data;
		std::vector<bool> exists_array;
	public:
		collocation_grid (int i_derivs, int i_rows, int i_cols) {rows = i_rows; cols = i_cols; derivs = i_derivs; data.resize (i_rows * i_cols * i_derivs); exists_array.resize (i_rows * i_cols * i_derivs, false);}
		virtual ~collocation_grid () {};
		inline double& operator () (int deriv, int row, int col) {index (deriv, row, col);}
		inline double& index (int deriv, int row, int col) {data [deriv * rows * cols + row * cols + col];}
		inline double& exists (int deriv, int row, int col) {exists_array [deriv * rows * cols + row * cols + col];}
	};
	
	class cheb_grid : public collocation_grid
	{
	private:
		double pioN;
		inline double recursion (int d, int m, int k) {
			if (exists (d, m, k)) {
				return index (d, m, k);
			} else if ((d == 2 && (m == 0 || m == 1) || (d == 1 && m == 0))) {
				return 0.0;
			} else if ((d == 1 && m == 1) || (d == 0 && m == 0)) {
				return 1.0;
			} else if (d == 0 && m == 1) {
				return pioN * k;
			} else if (d == 0) {
				return 2.0 * pioN * k * recursion (0, m - 1, k) + recursion (0, m - 2, k);
			} else if (d == 1) {
				return 2.0 * recursion (0, m - 1, k) + 2.0 * pioN * k * recursion (1, m - 1, k) - recursion (1, m - 2, k);
			} else if (d == 2) {
				return 4.0 * recursion (1, m - 1, k) + 2.0 * pioN * k * recursion (2, m - 1, k) - recursion (2, m - 2, k);
			} else {
				throw 0;
			}
		}
	public:
		cheb_grid (int i_M, int i_N) : collocation_grid (3, i_M, i_N) {
			int d, m, k;
			pioN = std::acos (-1.0) / i_N;
			for (d = 0; d < 3; ++d) {
				for (m = 0; m < i_M; ++m) {
					for (k = 0; k < i_N; ++k) {
						index (d, m, k) = recursion (d, m, k);
						exists (d, m, k) = true;
					}
				}
			}
		}
		virtual ~cheb_grid () {};
	};
} /* collocation */

#endif /* end of include guard: COLLOCATION_HPP_HLJUSAEZ */
