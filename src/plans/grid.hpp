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
#include <cmath>
#include "logger/logger.hpp"

namespace plans
{
	/*!**********************************************************************
	 * \brief A class containing the basic dimension information
	 ************************************************************************/
	class axis
	{
	private:
		int n; //!< The integer number of meshpoints along the axis
		int excess_0; //!< The integer number of excess meshpoints around index 0
		int excess_n; //!< The integer number of excess meshpoints around index n - 1
		double position_0; //!< The double position value at index excess_0
		double position_n; //!< The double position value at index n - 1 - excess_n
		
	public:
		/*!**********************************************************************
		 * \param i_n The integer number of meshpoints along the axis
		 * \param i_position_0 The double position value at index i_excess_0
		 * \param i_position_n The double position value at index n - 1 - i_excess_n
		 * \param i_excess_0 The integer number of excess meshpoints around index 0
		 * \param i_excess_n The integer number of excess meshpoints around index n - 1
		 ************************************************************************/
		axis (int i_n = 0, double i_position_0 = 0.0, double i_position_n = 0.0, int i_excess_0 = 0, int i_excess_n = 0) :
		n (i_n), excess_0 (i_excess_0), excess_n (i_excess_n), position_0 (i_position_0), position_n (i_position_n) {}
		
		/*!**********************************************************************
		 * \brief Get the total number of gridpoints along the axis
		 ************************************************************************/
		int get_n () {
			return n;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points near index 0
		 ************************************************************************/
		int get_excess_0 () {
			return excess_0;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points near index n - 1
		 ************************************************************************/
		int get_excess_n () {
			return excess_n;
		}
		
		/*!**********************************************************************
		 * \brief Get the position at index excess_0
		 ************************************************************************/
		double get_position_0 () {
			return position_0;
		}
		
		/*!**********************************************************************
		 * \brief Get the position at index n - 1 - excess_n
		 ************************************************************************/
		double get_position_n () {
			return position_n;
		}
	};
	
	/*!*******************************************************************
	 * \brief A class containing a collocation grid
	 * 
	 * This will need to be instantiated for a collocation method. Ideally, the structure contains a grid of polynomial values at a number of collocation points which can be called for little temporal expense. Unfortunately, creating such a grid is an expensive procedure, so the creation is delayed until the get_data method is called.
	 *********************************************************************/
	template <class datatype>
	class grid
	{
	protected:
		int n; //!< The number of meshpoints along the axis
		int ld; //!< The number of meshpoints in the data along the axis (accounts for buffer points)
		int excess_0; //!< The number of excess meshpoints near index 0
		int excess_n; //!< The number of excess meshpoints near index n - 1
		datatype position_0; //!< The position at index excess_0
		datatype position_n; //!< The position at index n - 1 - excess_n
		int derivs; //!< The integer number of derivatives deep the collocation grid runs
		
		std::vector <datatype> positions; //!< A vector containing the positions along the axis
		bool calculated_matrix; //!< A boolean describing whether the collocation matrix has been calculated
	
	private:
		std::vector<std::vector<datatype> > data; //!< A double vector containing the vectors of collocation data

	public:
		/*!*******************************************************************
		 * \param i_axis_ptr A pointer to an axis object
		 * \param i_derivs The integer number of derivatives in the grid
		 * \param i_ld
		 *********************************************************************/
		grid (axis *i_axis_ptr, int i_derivs, int i_ld = 0) :
		n (i_axis_ptr->get_n ()),
		ld (i_ld),
		excess_0 (i_axis_ptr->get_excess_0 ()),
		excess_n (i_axis_ptr->get_excess_n ()),
		position_0 (i_axis_ptr->get_position_0 ()),
		position_n (i_axis_ptr->get_position_n ()),
		derivs (i_derivs) {
			TRACE ("Instantiating...");
			if (ld == 0) {
				ld = n;
			}
			calculated_matrix = false;

			if (position_0 == position_n || n == 0) {
				ERROR ("Grid has no physical extent.");
				throw 0;
			}
			positions.resize (n + 1);

			TRACE ("Instantiated...");
		}
	
		virtual ~grid () {}
		
		/*!**********************************************************************
		 * \brief Return a reference to the position along the grid
		 * 
		 * \param index The integer index of interest
		 * 
		 * \return A reference to a position along the grid
		 ************************************************************************/
		const datatype& operator[] (int index) {
			return positions [index];
		}
		
		/*!**********************************************************************
		 * \brief Get the number of meshpoints along the axis
		 ************************************************************************/
		int get_n () {
			return n;
		}
		
		/*!**********************************************************************
		 * \brief Get the total array length along the axis (allowing for buffer points)
		 ************************************************************************/
		int get_ld () {
			return ld;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points by index 0
		 ************************************************************************/
		int get_excess_0 () {
			return excess_0;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points by index n - 1
		 ************************************************************************/
		int get_excess_n () {
			return excess_n;
		}
	
		/*!*******************************************************************
		 * \brief Get the data array for a given derivative level
		 * 
		 * \param deriv An integer derivative level (0 = value, 1 = first derivative, ...)
		 * 
		 * 
		 * This method is intended for use when a subroutine requires the matrix containing the collocation elements. If the matrix has not yet been calculated, this method will do so.
		 * 
		 * \return A pointer to the first element of the double data array
		 *********************************************************************/
		inline datatype* get_data (int deriv) {
			if (!calculated_matrix) {
				calculate_matrix ();
			}
			return &(data [deriv] [0]);
		}
		
	protected:
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
	
		/*!**********************************************************************
		 * \brief Calculate the collocation grid
		 * 
		 * This must be overwritten in subclasses that implement a grid.
		 ************************************************************************/
		virtual void _calculate_matrix () = 0;
		
	private:
		/*!**********************************************************************
		 * \brief Instantiate the collocation grid
		 * 
		 * Instantiate the collocation grid to be used by a collocation method. This method depends on the _calculate_matrix () operation.
		 ************************************************************************/
		virtual void calculate_matrix () {
			data.resize (derivs);
			for (int i = 0; i < derivs; ++i) {
				data [i].resize (ld * ld);
			}
			_calculate_matrix ();
			calculated_matrix = true;
		}
	};
	
#ifndef _VCOS
	enum mode {
		mode_flag = 0x10
	};

	namespace vertical
	{
		/*!*******************************************************************
		 * \brief A collocation grid for Chebyshev polynomials
		 * 
		 * This collocation grid stores the N collocation points for up to the Mth order Chebyshev polynomial and their first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public plans::grid <datatype>
		{
		private:
			using plans::grid <datatype>::n;
			using plans::grid <datatype>::ld;
			using plans::grid <datatype>::excess_0;
			using plans::grid <datatype>::excess_n;
			using plans::grid <datatype>::position_0;
			using plans::grid <datatype>::position_n;
			using plans::grid <datatype>::derivs;
			using plans::grid <datatype>::positions;
	
			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
			
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
	
			virtual ~grid () {};
			
		protected:
			/*!**********************************************************************
			 * \copydoc plans::grid::_calculate_matrix()
			 ************************************************************************/
			void _calculate_matrix ();
			
		private:
			/*!*******************************************************************
			* \brief A check to see whether an index exists, for recursion
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
			 * \param d The integer deriv to be indexed
			 * \param m The integer row to be indexed
			 * \param k The integer column to be indexed
			 * 
			 * Uses known relationships of Chebyshev polynomials to calculate the value of the polynomial and its derivatives at an arbitrary collocation grid point.
			 * 
			 * \return The datatype value of the index
			 *********************************************************************/
			datatype recursion (int d, int m, int k);
		};
	} /* vertical */
#else
	enum mode {
		mode_flag = 0x20
	};
	
	namespace vertical
	{
		/*!*******************************************************************
		 * \brief A collocation grid for Chebyshev polynomials
		 *
		 * This collocation grid stores the N collocation points for up to the Mth order cosine modes and their first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public plans::grid <datatype>
		{
		private:
			using plans::grid <datatype>::n;
			using plans::grid <datatype>::ld;
			using plans::grid <datatype>::excess_0;
			using plans::grid <datatype>::excess_n;
			using plans::grid <datatype>::position_0;
			using plans::grid <datatype>::position_n;
			using plans::grid <datatype>::derivs;
			using plans::grid <datatype>::positions;

			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations

		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);

			virtual ~grid () {};

		protected:
			/*!**********************************************************************
			 * \copydoc plans::grid::_calculate_matrix()
			 ************************************************************************/
			void _calculate_matrix ();
		};
	} /* vertical */
#endif
	
	namespace horizontal
	{
		/*!*******************************************************************
		 * \brief A collocation grid for Chebyshev polynomials
		 * 
		 * This collocation grid stores the N collocation points for up to the Mth order fourier mode and their first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public plans::grid <datatype>
		{
		private:
			using plans::grid <datatype>::n;
			using plans::grid <datatype>::ld;
			using plans::grid <datatype>::excess_0;
			using plans::grid <datatype>::excess_n;
			using plans::grid <datatype>::position_0;
			using plans::grid <datatype>::position_n;
			using plans::grid <datatype>::derivs;
			using plans::grid <datatype>::positions;
	
			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
			
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
				
			virtual ~grid () {};
			
		protected:
			/*!**********************************************************************
			 * \copydoc plans::grid::_calculate_matrix()
			 ************************************************************************/
			void _calculate_matrix ();
		};
	} /* fourier */
} /* plans */

#endif /* end of include guard: COLLOCATION_HPP_3FRTUK5Z */
