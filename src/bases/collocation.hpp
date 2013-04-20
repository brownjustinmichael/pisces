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
} /* bases */

#endif /* end of include guard: COLLOCATION_HPP_3FRTUK5Z */
