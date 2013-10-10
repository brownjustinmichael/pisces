/*!***********************************************************************
 * \file bases/collocation.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_HPP_3FRTUK5Z
#define COLLOCATION_HPP_3FRTUK5Z

#include <vector>
#include <memory>
#include "../config.hpp"

namespace bases
{
	/*!*******************************************************************
	 * \brief A class containing a collocation grid
	 * 
	 * This will need to be instantiated for a collocation method. 
	 * Ideally, the structure contains a grid of polynomial values at a 
	 * number of collocation points which can be called for little 
	 * temporal expense.
	 *********************************************************************/
	template <class datatype>
	class collocation_grid
	{
	public:
		/*!*******************************************************************
		 * \param i_derivs The integer number of derivatives in the grid
		 * \param i_rows The integer number of rows in the grid
		 * \param i_cols The integer number of columns in the grid
		 *********************************************************************/
		collocation_grid (int i_derivs, int i_rows, int i_cols) {
			rows = i_rows;
			cols = i_cols;
			derivs = i_derivs;

			TRACE ("Instantiating...")

			data.resize (derivs);

			for (int i = 0; i < i_derivs; ++i) {
				data [i].resize (i_rows * i_cols);
			}

			TRACE ("Instantiated...")
		}
	
		virtual ~collocation_grid () {}
	
		/*!*******************************************************************
		 * \brief An indexing operation into the grid, for convenience
		 * 
		 * \param deriv The integer deriv to be indexed
		 * \param row The integer row to be indexed
		 * \param col The integer column to be indexed
		 * 
		 * \return The double value at the index
		 *********************************************************************/
		inline datatype& index (int deriv, int row, int col) {
			return data [deriv] [row * cols + col];
		}
	
		/*!*******************************************************************
		 * \brief Get the data array for a given derivative level
		 * 
		 * This method is intended for use when a subroutine requires the 
		 * matrix containing the collocation elements.
		 * 
		 * \param deriv An integer derivative level (0 = value, 1 = first derivative, ...)
		 * 
		 * \return A pointer to the first element of the double data array
		 *********************************************************************/
		inline datatype* get_data (int deriv) {
			return &(data [deriv] [0]);
		}
	
	protected:
		int rows; //!< The integer number of rows in the grid
		int cols; //!< The integer number of columns in the grid
		int derivs; //!< The integer number of derivatives deep the collocation grid runs
	
	private:
		std::vector<std::vector<datatype> > data; //!< A double vector containing the collocation points
	};
	
	namespace chebyshev
	{
	   /*!*******************************************************************
	    * \brief A collocation grid for Chebyshev polynomials
	    * 
	    * This collocation grid stores the N collocation points for up to the
	    * Mth order Chebyshev polynomial and its first and second derivatives
	    *********************************************************************/
	   template <class datatype>
	   class grid : public bases::collocation_grid <datatype>
	   {
	   public:
	   	/*!*******************************************************************
	   	 * \param i_M The integer max order of Chebyshev polynomial
	   	 * \param i_N The integer number of collocation points
	   	 * \param i_scale A datatype by which the grid should be scaled
	   	 * \param i_width The datatype width of the collocation region
	   	 *********************************************************************/
	   	grid (int i_M, int i_N, datatype i_scale = 1.0, datatype i_width = 2.0);
	
	   	virtual ~grid () {};
	
	   private:
	   	using bases::collocation_grid <datatype>::derivs;
	   	using bases::collocation_grid <datatype>::rows;
	   	using bases::collocation_grid <datatype>::cols;
	
	   	datatype scale; //!< A datatype by which the collocation grid should be scaled
	   	datatype width; //!< The datatype width of the collocation region
	   	std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
	   	datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
	
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
	   	 * \return The datatype value of the index
	   	 *********************************************************************/
	   	datatype recursion (int d, int m, int k);
	   };
	} /* chebyshev */
	
	namespace fourier
	{
	   /*!*******************************************************************
	    * \brief A collocation grid for Chebyshev polynomials
	    * 
	    * This collocation grid stores the N collocation points for up to the
	    * Mth order Chebyshev polynomial and its first and second derivatives
	    *********************************************************************/
	   template <class datatype>
	   class grid : public bases::collocation_grid <datatype>
	   {
	   public:
	   	/*!*******************************************************************
	   	 * \param i_M The integer max order of Chebyshev polynomial
	   	 * \param i_N The integer number of collocation points
	   	 * \param i_scale A datatype by which the grid should be scaled
	   	 * \param i_width The datatype width of the collocation region
	   	 *********************************************************************/
	   	grid (int i_M, int i_N, datatype i_scale = 1.0, datatype i_width = 2.0);
	
	   	virtual ~grid () {};
	
	   private:
	   	using bases::collocation_grid <datatype>::derivs;
	   	using bases::collocation_grid <datatype>::rows;
	   	using bases::collocation_grid <datatype>::cols;
	
	   	datatype scale; //!< A datatype by which the collocation grid should be scaled
	   	datatype width; //!< The datatype width of the collocation region
	   	std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
	   	datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
	   };
	} /* fourier */
} /* bases */

#endif /* end of include guard: COLLOCATION_HPP_3FRTUK5Z */
