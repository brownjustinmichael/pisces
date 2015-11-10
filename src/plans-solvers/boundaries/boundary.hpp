/*!**********************************************************************
 * \file boundary.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-07-14.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_013D6464
#define BOUNDARY_HPP_013D6464

#include "versions/version.hpp"
#include "plans/grids/grid.hpp"

/*!**********************************************************************
 * \namespace boundaries
 * 
 * \brief A namespace that contains boundary objects used to calculate matrices and right hand sides across boundaries
 ************************************************************************/
namespace boundaries
{
	/*
		TODO Make non ASCII documentation for this class
	*/

	/*!**********************************************************************
	 * \brief A class used to calculate cross-boundary matrices
	 * 
	 * The main members of this class are the calculate_rhs and calculate_matrix members, which act on the overlapping terms in a block matrix. It's important to understand the boundary setup:
	 * 
	 * Element 1: ---X-----X-----X-----B1--E
	 * Element 2:                  E---B2----X-----X-----X---
	 * 
	 * For stability, we have added a slight overlap from element to element. This ensures smooth derivatives at the boundaries. These overlapping poins do not contribute to the solve, but rather are forced to be a value interpolated from the solution of the adjacent element. A bit of terminology: the "boundary" is the redundant point that occurs in both elements B. The "overlap" region includes every B and E point. The "excess" is the non-redudant overlap E.
	 * 
	 * Put into matrix form, this looks like the following:
	 * 
	 * Element 1 Element 2
	 * 
	 * XXXXXXXXX          <- Inside Element 1
	 * BBBBBBBBB--------- <- Boundary of Element 1 (data pointers for bottom boundaries should point here)
	 * EEEEEEEEEEEEEEEEEE <- Excess region inside Element 1
	 * EEEEEEEEEEEEEEEEEE <- Excess region inside Element 2
	 * ---------BBBBBBBBB <- Boundary of Element 2 (data pointers for upper bounderies should pointe here)
	 *          XXXXXXXXX <- Inside Element 2
	 * 
	 * Note that the boundary point appears twice in the matrix. This is not an oversight; the two conditions to be met are 1) ensure that the boundary evolves according to the scheme and 2) ensure that B1 and B2 have the same value.
	 ************************************************************************/
	template <class datatype>
	class boundary
	{
	public:
		boundary () {}
		
		virtual ~boundary () {}

		/**
		 * @brief A factory from which a boundary instance can be generated
		 */
		class factory
		{
		public:			
			virtual ~factory () {}
			
			/**
			 * @brief Generate an instance of a boundary object associated with the factory
			 * 
			 * @param grids An array of grid objects needed to construct the boundary
			 * @param top Whether this is a top boundary
			 * 
			 * @return A shared pointer to the newly constructed boundary class
			 */
			virtual std::shared_ptr <boundary <datatype>> instance (grids::grid **grids, bool top) = 0;

			/**
			 * @brief Generate an instance of a boundary object associated with the factory
			 * 
			 * @param grid_n A reference to the horizontal grid
			 * @param grid_m A reference to the vertical grid
			 * @param top Whether this is a top boundary
			 * 
			 * @return A shared pointer to the newly constructed boundary class
			 */
			virtual std::shared_ptr <boundary <datatype>> instance (grids::grid &grid_n, grids::grid &grid_m, bool top) {
				grids::grid *grids [3] = {&grid_n, &grid_m};
				return instance (grids, top);
			}
		};
		
		/*!**********************************************************************
		 * \brief Get the version of the class
		 * 
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.1.0.0");
			return version;
		}
		
		/*!**********************************************************************
		 * \brief Get the total overlap region for this boundary
		 * 
		 * This member returns the number of matrix rows total that overlap between this element and the adjacent one
		 * 
		 * \return The integer overlap region
		 ************************************************************************/
		virtual int get_overlap () {
			return 0;
		}
		
		/*!**********************************************************************
		 * \brief Get the overlap region for this boundary outside the current element
		 * 
		 * This member returns how far this element penetrates into an adjacent one
		 * 
		 * \return The integer overlap region outside the element
		 ************************************************************************/
		virtual int get_ex_overlap () {
			return 0;
		}
		
		/*!**********************************************************************
		 * \brief Get the excess region for this boundary inside the current element excluding any redundancies
		 * 
		 * In general, there is one point that is contained in both elements, so this will usually return one less than the overlap region inside the boundary.
		 * 
		 * \return The integre excess region inside the element
		 ************************************************************************/
		virtual int get_excess () {
			return 0;
		}
		
		/*
			TODO Remove dimensionality from base class
		*/
		
		/*!**********************************************************************
		 * \brief If this is a communicating boundary, send data to the adjacent element, else do nothing
		 * 
		 * \param data_temp The datatype pointer to the data to send
		 * \param lda The leading dimension of the data to send
		 * \param n The depth of the data to send, if negative, send the overlap region
		 * 
		 * This sends data_temp to the adjacent element if it exists. This will automatically send the full width of the boundary and will send the first n rows of data_temp.
		 ************************************************************************/
		virtual void send (datatype *data_temp, int lda, int n = -1) {}
		
		/*!**********************************************************************
		 * \brief If this is a communicating boundary, receive data from the adjacent element, else do nothing
		 * 
		 * \param data_temp The datatype pointer to the data to receive
		 * \param lda The leading dimension of the data to receive
		 * \param n The depth of the data to receive, if negative, receive the overlap region
		 * \param alpha Scale data_temp by alpha before adding in the received data (0.0 for a clean copy)
		 * 
		 * This receives data_temp to the adjacent element if it exists. This will automatically receive the full width of the boundary and will receive the first n rows of data_temp.
		 ************************************************************************/
		virtual void receive (datatype *data_temp, int lda, int n = -1, datatype alpha = 1.0) {}
		
		/*!**********************************************************************
		 * \brief Calculate the right hand side of the matrix equation
		 * 
		 * \param data The pointer to the previous timestep data at the boundary
		 * \param data_temp The pointer to the right hand side data at the boundary
		 * \param m The leading dimension of data
		 * \param lda The leading dimension of data_temp
		 * \param flag A flag indicating whether this is an x_solve or z_solve
		 ************************************************************************/
		virtual void calculate_rhs (datatype *data, datatype *data_temp, int m, int lda, int flag) = 0;
		
		/*!**********************************************************************
		 * \brief Calculate the matrix for the matrix equation
		 * 
		 * \param timestep The timestep needed to construct the matrix
		 * \param default_matrix The matrix corresponding to the values of the data pointed at the internal boundary, assumed to have leading dimension m
		 * \param matrix_in The matrix containing the appropriate derivative information for the equation pointed at the internal boundary, assumed to have leading dimension m
		 * \param matrix_out The array to which the final matrix is written
		 * \param lda The leading dimension of matrix_out
		 * \param diverging Whether this matrix is associated with the incompressible corrector object
		 ************************************************************************/
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *matrix_out, int lda, bool diverging = false) = 0;
	};

	/**
	 * @brief Match the boundaries of a class if they are communicating boundaries
	 * @details This function was added as a shorthand to match any overlapping regions of elements quickly
	 * 
	 * @param m The number of datapoints in the vertical
	 * @param boundary_0 A shared pointer to the top boundary
	 * @param boundary_n A shared pointer to the bottom boundary
	 * @param data The pointer to the data 
	 * @param excess_0 The number of excess points past the top boundary
	 * @param excess_n The number of excess points past the bottom boundary
	 */
	template <class datatype>
	void boundary_match (int m, std::shared_ptr <boundaries::boundary <datatype>> &boundary_0, std::shared_ptr <boundaries::boundary <datatype>> &boundary_n, datatype *data, int excess_0 = 1, int excess_n = 1) {
	if (boundary_n) {
		boundary_n->send (data + m - 2 * excess_n - 1, m, excess_n);
	}
	if (boundary_0) {
		boundary_0->receive (data, m, excess_0, 0.0);
		boundary_0->send (data + excess_0, m, excess_0 + 1);
	}
	if (boundary_n) {
		boundary_n->receive (data + m - excess_n - 1, m, excess_n + 1, 0.0);
	}
}

} /* boundaries */

#endif /* end of include guard: BOUNDARY_HPP_013D6464 */

