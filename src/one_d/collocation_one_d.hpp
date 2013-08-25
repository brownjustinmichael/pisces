/*!***********************************************************************
 * \file collocation_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_HPP_HLJUSAEZ
#define COLLOCATION_HPP_HLJUSAEZ

#include <memory>
#include <vector>
#include <cmath>
#include "../bases/collocation.hpp"

namespace one_d
{
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
} /* one_d */


#endif /* end of include guard: COLLOCATION_HPP_HLJUSAEZ */
