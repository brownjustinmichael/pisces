/*!**********************************************************************
 * \file implemented_boundary.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-07-14.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_TWO_D_HPP_863EDB36
#define BOUNDARY_TWO_D_HPP_863EDB36

#include "mpi/messenger.hpp"
#include "linalg/utils.hpp"
#include "plans/grids/grid.hpp"
#include "plans-transforms/transform.hpp"
#include "plans-transforms/transformer.hpp"

#include "boundary.hpp"

namespace boundaries
{
	/*!**********************************************************************
	 * \brief A boundary that fixes the value of the data at the edge
	 ************************************************************************/
	template <class datatype>
	class fixed_boundary : public boundaries::boundary <datatype>
	{
	private:
		datatype value; //!< The value at which to fix the data
		bool top; //!< A boolean regarding whether this is the top or bottom of an element
		int ldn; //!< The horizontal array extent
		int m; //!< The number of data points in the vertical
		
	public:
		/*!**********************************************************************
		 * \param i_grid_n The horizontal grid object
		 * \param i_grid_m The vertical grid object
		 * \param i_value The value at which to fix the data at the edge
		 * \param i_top A boolean regarding whether this is the top or bottom of an element
		 ************************************************************************/
		fixed_boundary (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_value, bool i_top) : value (i_value * std::sqrt (i_grid_n.get_n ())), top (i_top) {
			ldn = i_grid_n.get_ld ();
			m = i_grid_m.get_n ();
		}
		
		virtual ~fixed_boundary () {}

		/**
		 * @copydoc boundary::factory
		 */
		class factory : public boundaries::boundary <datatype>::factory
		{
		private:
			datatype value; //!< The value of the boundary

		public:
			/**
			 * @param i_value The value of the boundary
			 */
			factory (datatype i_value) : value (i_value) {}

			virtual ~factory () {}
			
			/**
			 * @copydoc boundary::factory::instance
			 */
			virtual std::shared_ptr <boundary <datatype>> instance (grids::grid <datatype> **grids, bool top) {
				return std::shared_ptr <boundary <datatype>> (new fixed_boundary (*grids [0], *grids [1], value, top));
			}
		};
		
		/*!**********************************************************************
		 * \copydoc boundary::calculate_rhs
		 ************************************************************************/
		virtual void calculate_rhs (datatype *data, datatype *data_temp, int m, int lda, int flag) {
			// Because the equation will be solved in fourier space, set the first value to the desired value and the rest to 0
			data_temp [0] = value;
			for (int i = 1; i < ldn; ++i) {
				data_temp [i * lda] = 0.0;
			}
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::calculate_matrix
		 ************************************************************************/
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *matrix_out, int lda, bool diverging = false) {
			linalg::scale (m, 0.0, matrix_out, lda);
			linalg::add_scaled (m, 1.0, default_matrix, matrix_out, m, lda);
		}
	};
	
	/*!**********************************************************************
	 * \brief A boundary the fixes the derivative of the data at the edge
	 ************************************************************************/
	template <class datatype>
	class fixed_deriv_boundary : public boundaries::boundary <datatype>
	{
	private:
		datatype value; //!< The value at which to fix the derivative
		bool top; //!< A boolean regarding whether this is the top or bottom of an element
		int ldn; //!< The horizontal array extent
		int m; //!< The number of data points in the vertical
		datatype *deriv_matrix; //!< The matrix containing the derivative information
		
	public:
		/*!**********************************************************************
		 * \param i_grid_n The horizontal grid object
		 * \param i_grid_m The vertical grid object
		 * \param i_value The value at which to fix the derivative at the edge
		 * \param i_top A boolean regarding whether this is the top or bottom of an element
		 ************************************************************************/
		fixed_deriv_boundary (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype i_value, bool i_top) : value (i_value * std::sqrt (i_grid_n.get_n ())), top (i_top) {
			ldn = i_grid_n.get_ld ();
			m = i_grid_m.get_n ();
			
			// Get the derivative matrix from the grid object
			deriv_matrix = i_grid_m.get_data (1) + (top ? m - 1 : 0);
		}
		
		virtual ~fixed_deriv_boundary () {}

		/**
		 * @copydoc boundary::factory
		 */
		class factory : public boundaries::boundary <datatype>::factory
		{
		private:
			datatype value; //!< The value at which the derivative will be fixed

		public:
			/**
			 * @param i_value The value at which the derivative will be fixed
			 */
			factory (datatype i_value) : value (i_value) {}

			virtual ~factory () {}
			
			/**
			 * @copydoc boundary::factory::instance
			 */
			virtual std::shared_ptr <boundary <datatype>> instance (grids::grid <datatype> **grids, bool top) {
				return std::shared_ptr <boundary <datatype>> (new fixed_deriv_boundary (*grids [0], *grids [1], value, top));
			}
		};
		
		/*!**********************************************************************
		 * \copydoc boundary::calculate_rhs
		 ************************************************************************/
		virtual void calculate_rhs (datatype *data, datatype *data_temp, int m, int lda, int flag) {
			// If this is a z-directional solve, it's the same as the fixed_boundary case
			if (flag & z_solve) {
				data_temp [0] = value;
				for (int i = 1; i < ldn; ++i) {
					data_temp [i * lda] = 0.0;
				}
			}
			
			// If this is an x-directional solve, we have no knowledge of the z-derivative, so we fix the boundaries
			if (flag & x_solve) {
				for (int i = 0; i < ldn; ++i) {
					data_temp [i * lda] = data [i * m];
				}
			}
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::calculate_matrix
		 ************************************************************************/
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *matrix_out, int lda, bool diverging = false) {
			linalg::scale (m, 0.0, matrix_out, lda);
			linalg::add_scaled (m, 1.0, deriv_matrix, matrix_out, m, lda);
		}
	};
	
	/*!**********************************************************************
	 * \brief A boundary to communicate between adjacent elements
	 ************************************************************************/
	template <class datatype>
	class communicating_boundary : public boundaries::boundary <datatype>
	{
	private:
		int ldn; //!< The horizontal extent of the data array
		int m; //!< The number of data points in vertical in the data
		int n_boundary_in; //!< The number of internal excess points
		int n_boundary_out; //!< The number of external excess points
		
		mpi::messenger *messenger_ptr; //!< A pointer to the element messenger for communication
		int out_id; //!< The id of the adjacent element
		
		datatype alpha; //!< The alpha parameter (usually 0.5) indicating how we are winding across the boundary (0.0 would be adopt all properties from the external element)
		bool top; //!< A boolean regarding whether this is the top or bottom of an element
		
		std::vector <datatype> buffer; //!< A data buffer for convenience
		
	public:
		/*!**********************************************************************
		 * \param i_messenger_ptr A pointer to the element messenger for communication
		 * \param i_ldn The horizontal extent of the data array
		 * \param i_m The number of data points in vertical in the data
		 * \param i_n_boundary_in The number of internal excess points
		 * \param i_top A boolean regarding whether this is the top or bottom of an element
		 ************************************************************************/
		communicating_boundary (mpi::messenger *i_messenger_ptr, int i_ldn, int i_m, int i_n_boundary_in, bool i_top) : messenger_ptr (i_messenger_ptr) {
			/*
				TODO Take grid as input
			*/
			alpha = 0.5;
			ldn = i_ldn;
			m = i_m;
			n_boundary_in = i_n_boundary_in;
			top = i_top;
			
			// Communicate the excess information
			out_id = messenger_ptr->get_id () + (top ? 1 : -1);
			messenger_ptr->template send <int> (1, &i_n_boundary_in, out_id, 0);
			messenger_ptr->template recv <int> (1, &n_boundary_out, out_id, 0);
			
			buffer.resize ((std::max (n_boundary_out, n_boundary_in) + 1) * ldn);
		}
		
		virtual ~communicating_boundary () {}

		/**
		 * @copydoc boundary::factory
		 */
		class factory : public boundaries::boundary <datatype>::factory
		{
		private:
			mpi::messenger *messenger_ptr; //!< A pointer to the mpi messenger object

		public:
			/**
			 * @param i_messenger_ptr A pointer to the mpi messenger object
			 */
			factory (mpi::messenger *i_messenger_ptr) : messenger_ptr (i_messenger_ptr) {}

			virtual ~factory () {}
			
			/**
			 * @copydoc boundary::factory::instance
			 */
			virtual std::shared_ptr <boundary <datatype>> instance (grids::grid <datatype> **grids, bool top) {
				return std::shared_ptr <boundary <datatype>> (new communicating_boundary (messenger_ptr, grids [0]->get_ld (), grids [1]->get_n (), top ? grids [1]->get_excess_n () : grids [1]->get_excess_0 (), top));
			}
		};
		
		/*!**********************************************************************
		 * \copydoc boundary::get_overlap
		 ************************************************************************/
		virtual int get_overlap () {
			return 2 + n_boundary_in + n_boundary_out;
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::get_ex_overlap
		 ************************************************************************/
		virtual int get_ex_overlap () {
			return 1 + n_boundary_out;
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::get_excess
		 ************************************************************************/
		virtual int get_excess () {
			return n_boundary_in;
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::send
		 ************************************************************************/
		virtual void send (datatype *data_temp, int lda, int n = -1) {
			for (int i = 0; i < ldn; ++i) {
				for (int j = 0; j < (n < 0 ? n_boundary_out + 1 : n); ++j) {
					buffer [j * ldn + i] = data_temp [i * lda + j];
				}
			}

			messenger_ptr->send ((n_boundary_out + 1) * ldn, &buffer [0], out_id, 0);
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::receive
		 ************************************************************************/
		virtual void receive (datatype *data_temp, int lda, int n = -1, datatype alpha = 1.0) {
			messenger_ptr->recv ((n_boundary_in + 1) * ldn, &buffer [0], out_id, 0);
			
			for (int i = 0; i < ldn; ++i) {
				for (int j = 0; j < (n < 0 ? n_boundary_out + 1 : n); ++j) {
					// DEBUG ("PREV: " << data_temp [i * lda + j]);
					data_temp [i * lda + j] = alpha * data_temp [i * lda + j] + buffer [j * ldn + i];
				}
			}
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::calculate_rhs
		 ************************************************************************/
		virtual void calculate_rhs (datatype *data, datatype *data_temp, int m, int lda, int flag) {
			// Zero everything but the internal boundary row
			linalg::matrix_scale (1 + n_boundary_out + n_boundary_in, ldn, 0.0, data_temp + (top ? 1 : -1 - n_boundary_out - n_boundary_in), lda);
			// Scale the internal boundary by alpha
			linalg::scale (ldn, alpha, data_temp, lda);
			// Add the original data to the internal boundary row, scaled by alpha (to ensure that the boundary averages out between the two elements)
			if (data) {
				linalg::add_scaled (ldn, alpha, data, data_temp, m, lda);
			}
			// Copy the internal boundary row to the external boundary row
			linalg::copy (ldn, data_temp, data_temp + (top ? 1 : -1) * (1 + n_boundary_out + n_boundary_in), lda, lda);
		}
		
		/*!**********************************************************************
		 * \copydoc boundary::calculate_matrix
		 ************************************************************************/
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *matrix_out, int lda, bool diverging = false) {
			// Zero everything but the internal boundary row
			linalg::matrix_scale (1 + n_boundary_out + n_boundary_in, m, 0.0, matrix_out + (top ? 1 : -1 - n_boundary_out - n_boundary_in), lda);
			// Setting the external boundary matrix row with the row from the internal boundary row
			linalg::matrix_add_scaled (1, m, alpha, matrix_in, matrix_out + (top ? 1 : -1) * (1 + n_boundary_out + n_boundary_in), m, lda);
			// Scale the new derivative information in the matrix by the timestep
			linalg::matrix_scale (1 + n_boundary_out + n_boundary_in + 1, m, timestep, matrix_out + (top ? 0 : -1 - n_boundary_out - n_boundary_in), lda);
			
			// For the internal excess region and the internal boundary, add the final data value
			linalg::matrix_add_scaled (1 + n_boundary_in, m, 1.0, default_matrix + (top ? 0 : -n_boundary_in), matrix_out + (top ? 0 : -n_boundary_in), m, lda);
			// For the external excess region, subtract the final data value
			// This ensures that the values of the excess region will exactly match the corresponding values in the adjacent element they're attempting to represent
			// A note: this is not strictly accurate since the points may not be at exactly the same position, but this is much faster and simpler than interpolation and is not far enough from the real value to justify the more complicated method
			linalg::matrix_add_scaled (n_boundary_out, m, -1.0, default_matrix + (top ? -1 : 1), matrix_out + (top ? 1 + n_boundary_in : -n_boundary_in - n_boundary_out), m, lda);
		}
	};
} /* boundaries */

#endif /* end of include guard: BOUNDARY_TWO_D_HPP_863EDB36 */
