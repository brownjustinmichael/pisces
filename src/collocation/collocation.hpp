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
		collocation_grid (int i_derivs, int i_rows, int i_cols);
		virtual ~collocation_grid () {};
		inline double& operator () (int deriv, int row, int col) {index (deriv, row, col);}
		double &index (int deriv, int row, int col);
		std::vector<bool>::reference exists (int deriv, int row, int col) {
			return exists_array [deriv * rows * cols + row * cols + col];
		}
	};
	
	class cheb_grid : public collocation_grid
	{
	private:
		double pioN;
		double recursion (int d, int m, int k);
	public:
		cheb_grid (int i_M, int i_N);
		virtual ~cheb_grid () {};
	};
} /* collocation */

#endif /* end of include guard: COLLOCATION_HPP_HLJUSAEZ */
