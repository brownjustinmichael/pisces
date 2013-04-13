/*!***********************************************************************
 * \file collocation.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_HPP_HLJUSAEZ
#define COLLOCATION_HPP_HLJUSAEZ

#include <vector>
#include <cmath>

namespace collocation
{
	/*!*******************************************************************
	 * \brief An abstract class containing a collocation grid
	 * 
	 * This will need to be instantiated for a collocation method. 
	 * Ideally, the structure contains a grid of polynomial values at a 
	 * number of collocation points which can be called for little 
	 * temporal expense.
	 *********************************************************************/
	class collocation_grid
	{
	public:
		/*!*******************************************************************
		 * \param i_derivs The integer number of derivatives in the grid
		 * \param i_rows The integer number of rows in the grid
		 * \param i_cols The integer number of columns in the grid
		 *********************************************************************/
		collocation_grid (int i_derivs, int i_rows, int i_cols);
		virtual ~collocation_grid () {};
		/*!*******************************************************************
		 * \brief An indexing operation into the grid, for convenience
		 * 
		 * \param deriv The integer deriv to be indexed
		 * \param row The integer row to be indexed
		 * \param col The integer column to be indexed
		 * 
		 * \return The double value at the index
		 *********************************************************************/
		inline double &index (int deriv, int row, int col) {
			return data [deriv] [row * cols + col];
		}
		inline double *get_data (int deriv) {
			return &(data [deriv] [0]);
		}
		
	protected:
		int rows; //!< The integer number of rows in the grid
		int cols; //!< The integer number of columns in the grid
		int derivs; //!< The integer number of derivatives deep the collocation grid runs
		
	private:
		std::vector<std::vector<double>> data; //!< A double vector containing the collocation points
	};
	
	/*!*******************************************************************
	 * \brief A collocation grid for Chebyshev polynomials
	 * 
	 * This collocation grid stores the N collocation points for up to the
	 * Mth order Chebyshev polynomial and its first and second derivatives
	 *********************************************************************/
	class chebyshev_grid : public collocation_grid
	{
	public:
		/*!*******************************************************************
		 * \param i_M The integer max order of Chebyshev polynomial
		 * \param i_N The integer number of collocation points
		 *********************************************************************/
		chebyshev_grid (int i_M, int i_N, double i_scale = 1.0);
		virtual ~chebyshev_grid () {};
		
	private:
		double scale;
		std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
		double pioN; //!< The double 3.14159.../N, for use in calculations
		/*!*******************************************************************
	 	* \brief A check to see whether an element exists, for recursion
	 	* 
		* \param d The integer deriv to be indexed
		* \param m The integer row to be indexed
		* \param k The integer column to be indexed
		* 
		* \return An std::vector<bool>::reference of whether the element exists
	 	*********************************************************************/
		std::vector<bool>::reference exists (int d, int m, int k) {
			return exists_array [d + m * derivs + k * derivs * rows];
		}
		/*!*******************************************************************
		 * \brief A recursive method to calculate the full collocation grid
		 * 
		 * Uses known relationships of Chebyshev polynomials to calculate the 
		 * value of the polynomial and its derivatives at an arbitrary 
		 * collocation grid point.
		 * 
		 * \param d The integer deriv to be indexed
		 * \param m The integer row to be indexed
		 * \param k The integer column to be indexed
		 * 
		 * \return The double value of the index
		 *********************************************************************/
		double recursion (int d, int m, int k);
		
	};
} /* collocation */

#endif /* end of include guard: COLLOCATION_HPP_HLJUSAEZ */
