/*!***********************************************************************
 * \file bases/grid.hpp
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
	/*!**********************************************************************
	 * \brief A class containing the basic dimension information
	 ************************************************************************/
	class axis
	{
	public:
		/*!**********************************************************************
		 * \param i_n The integer number of meshpoints along the axis
		 * \param i_position_0 The double position value at index i_excess_0
		 * \param i_position_n The double position value at index n - 1 - i_excess_n
		 * \param i_excess_0 The integer number of excess meshpoints around index 0
		 * \param i_excess_n The integer number of excess meshpoints around index n - 1
		 ************************************************************************/
		axis (int i_n, double i_position_0, double i_position_n, int i_excess_0 = 0, int i_excess_n = 0) :
		n (i_n), excess_0 (i_excess_0), excess_n (i_excess_n), position_0 (i_position_0), position_n (i_position_n) {}
		
		int n; //!< The integer number of meshpoints along the axis
		int excess_0; //!< The double position value at index excess_0
		int excess_n; //!< The double position value at index n - 1 - i_excess_n
		double position_0; //!< The integer number of excess meshpoints around index 0
		double position_n; //!< The integer number of excess meshpoints around index n - 1
	};
	
	/*!*******************************************************************
	 * \brief A class containing a collocation grid
	 * 
	 * This will need to be instantiated for a collocation method. 
	 * Ideally, the structure contains a grid of polynomial values at a 
	 * number of collocation points which can be called for little 
	 * temporal expense.
	 *********************************************************************/
	template <class datatype>
	class grid
	{
	public:
		/*!*******************************************************************
		 * \param i_axis_ptr A pointer to an axis object
		 * \param i_derivs The integer number of derivatives in the grid
		 *********************************************************************/
		grid (axis *i_axis_ptr, int i_derivs, int i_ld = 0) :
		n (i_axis_ptr->n),
		ld (i_ld),
		excess_0 (i_axis_ptr->excess_0),
		excess_n (i_axis_ptr->excess_n),
		position_0 (i_axis_ptr->position_0),
		position_n (i_axis_ptr->position_n),
		derivs (i_derivs) {
			TRACE ("Instantiating...");
			if (ld == 0) {
				ld = n;
			}

			data.resize (derivs);
			positions.resize (n + 1);

			for (int i = 0; i < i_derivs; ++i) {
				data [i].resize (ld * ld);
			}

			TRACE ("Instantiated...");
		}
		
		/*
			TODO If printing grids can be done logically, it would be nice to put cell in here, too.
		*/
	
		virtual ~grid () {
			// printf ("Destroying base grid\n");
		}
		
		/*!**********************************************************************
		 * \brief Return a reference to the position along the grid
		 * 
		 * \param index The integer index of interest
		 * 
		 * \return A reference to a position along the grid
		 ************************************************************************/
		datatype& operator[] (int index) {
			return positions [index];
		}
		
		/*!*******************************************************************
		 * \brief An indexing operation into the grid, for convenience
		 * 
		 * \param deriv The integer deriv to be indexed
		 * \param row The integer row to be indexed
		 * \param col The integer column to be indexed
		 * 
		 * \return The value at the index
		 *********************************************************************/
		inline datatype& index (int deriv, int row, int col) {
			return data [deriv] [row * ld + col];
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

		int n, ld; //!< The number of meshpoints along the axis
		int excess_0; //!< The number of excess meshpoints near index 0
		int excess_n; //!< The number of excess meshpoints near index n - 1
		datatype position_0; //!< The position at index excess_0
		datatype position_n; //!< The position at index n - 1 - excess_n
	
	protected:
		std::vector <datatype> positions; //!< A vector containing the positions along the axis
		int derivs; //!< The integer number of derivatives deep the collocation grid runs
	
	private:
		std::vector<std::vector<datatype>> data; //!< A double vector containing the vectors of collocation data
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
		class grid : public bases::grid <datatype>
		{
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
	
			virtual ~grid () {};
			
		private:
			using bases::grid <datatype>::n;
			using bases::grid <datatype>::ld;
			using bases::grid <datatype>::excess_0;
			using bases::grid <datatype>::excess_n;
			using bases::grid <datatype>::position_0;
			using bases::grid <datatype>::position_n;
			using bases::grid <datatype>::derivs;
			using bases::grid <datatype>::positions;
	
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
				return exists_array [d + m * derivs + k * derivs * n];
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
	
	namespace cosine
	{
		/*!*******************************************************************
		 * \brief A collocation grid for Chebyshev polynomials
		 * 
		 * This collocation grid stores the N collocation points for up to the
		 * Mth order Chebyshev polynomial and its first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public bases::grid <datatype>
		{
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
	
			virtual ~grid () {};
			
		private:
			using bases::grid <datatype>::n;
			using bases::grid <datatype>::ld;
			using bases::grid <datatype>::excess_0;
			using bases::grid <datatype>::excess_n;
			using bases::grid <datatype>::position_0;
			using bases::grid <datatype>::position_n;
			using bases::grid <datatype>::derivs;
			using bases::grid <datatype>::positions;
	
			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
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
		class grid : public bases::grid <datatype>
		{
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
	
			virtual ~grid () {};
			
		private:
			using bases::grid <datatype>::n;
			using bases::grid <datatype>::ld;
			using bases::grid <datatype>::excess_0;
			using bases::grid <datatype>::excess_n;
			using bases::grid <datatype>::position_0;
			using bases::grid <datatype>::position_n;
			using bases::grid <datatype>::derivs;
			using bases::grid <datatype>::positions;
	
			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
		};
	} /* fourier */
} /* bases */

#endif /* end of include guard: COLLOCATION_HPP_3FRTUK5Z */
