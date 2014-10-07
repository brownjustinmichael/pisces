/*!**********************************************************************
 * \file boundary_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-07-14.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_TWO_D_HPP_863EDB36
#define BOUNDARY_TWO_D_HPP_863EDB36

#include "mpi/messenger.hpp"
#include "linalg/utils.hpp"
#include "linalg/interpolate.hpp"
#include "plans/grid.hpp"

#include "boundary.hpp"
#include "solver.hpp"

namespace plans
{
	template <class datatype>
	class fixed_boundary : public plans::boundary <datatype>
	{
	private:
		datatype value;
		bool top;
		int ldn, m;
		
	public:
		fixed_boundary (plans::grid <datatype> *i_grid_n, plans::grid <datatype> *i_grid_m, datatype i_value, bool i_top) : value (i_value * std::sqrt (i_grid_n->get_n ())), top (i_top) {
			DEBUG ("Fixed boundary initialized");
			ldn = i_grid_n->get_ld ();
			m = i_grid_m->get_n ();
		}
		
		virtual ~fixed_boundary () {}
		
		virtual void calculate_rhs (datatype *data, datatype *interpolate_original, datatype *interpolate_data, datatype *data_temp, int m, int lda, int flag) {
			data_temp [0] = value;
			for (int i = 1; i < ldn; ++i) {
				data_temp [i * lda] = 0.0;
			}
		}
		
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *interpolate_matrix, datatype *matrix_out, int lda, bool diverging = false) {
			linalg::scale (m, 0.0, matrix_out, lda);
			linalg::add_scaled (m, 1.0, default_matrix, matrix_out, m, lda);
		}
	};
	
	template <class datatype>
	class fixed_deriv_boundary : public plans::boundary <datatype>
	{
	private:
		datatype value;
		bool top;
		int ldn, m;
		datatype *deriv_matrix;
		
	public:
		fixed_deriv_boundary (plans::grid <datatype> *i_grid_n, plans::grid <datatype> *i_grid_m, datatype i_value, bool i_top) : value (i_value * std::sqrt (i_grid_n->get_n ())), top (i_top) {
			DEBUG ("Fixed boundary initialized");
			ldn = i_grid_n->get_ld ();
			m = i_grid_m->get_n ();
			deriv_matrix = i_grid_m->get_data (1) + (top ? m - 1 : 0);
			DEBUG ("Done");
		}
		
		virtual ~fixed_deriv_boundary () {}
		
		virtual void calculate_rhs (datatype *data, datatype *interpolate_original, datatype *interpolate_data, datatype *data_temp, int m, int lda, int flag) {
			DEBUG ("Call " << data_temp);
			if (flag & z_solve) {
				data_temp [0] = value;
				for (int i = 1; i < ldn; ++i) {
					data_temp [i * lda] = 0.0;
				}
			}
			if (flag & x_solve) {
				for (int i = 0; i < ldn; ++i) {
					data_temp [i * lda] = data [i * m];
				}
			}
		}
		
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *interpolate_matrix, datatype *matrix_out, int lda, bool diverging = false) {
			linalg::scale (m, 0.0, matrix_out, lda);
			linalg::add_scaled (m, 1.0, deriv_matrix, matrix_out, m, lda);
		}
	};
	
	template <class datatype>
	class communicating_boundary : public plans::boundary <datatype>
	{
	private:
		datatype alpha;
		mpi::messenger *messenger_ptr;
		bool top;
		int n_boundary_in;
		int n_boundary_out;
		int out_id;
		const datatype *positions;
		datatype *boundary_positions;
		std::vector <datatype> boundary_positions_vec;
		std::vector <datatype> buffer;
		int ldn, m;
		
	public:
		communicating_boundary (mpi::messenger *i_messenger_ptr, int i_ldn, int i_m, int i_n_boundary_in, const datatype *i_positions, int index, bool i_top) : messenger_ptr (i_messenger_ptr) {
			/*
				TODO Take grid as input
			*/
			DEBUG ("Boundary Initialized");
			alpha = 0.5;
			positions = i_positions;
			n_boundary_in = i_n_boundary_in;
			top = i_top;
			out_id = messenger_ptr->get_id () + (top ? 1 : -1);
			messenger_ptr->template send <int> (1, &i_n_boundary_in, out_id, 0);
			messenger_ptr->template recv <int> (1, &n_boundary_out, out_id, 0);
			boundary_positions_vec.resize (n_boundary_out);
			boundary_positions = &boundary_positions_vec [0];
			messenger_ptr->template send <datatype> (i_n_boundary_in, i_positions + index, out_id, 1);
			messenger_ptr->template recv <datatype> (n_boundary_out, boundary_positions, out_id, 1);
			ldn = i_ldn;
			m = i_m;
		}
		
		virtual ~communicating_boundary () {}
		
		virtual int get_overlap () {
			return 2 + n_boundary_in + n_boundary_out;
		}
		
		virtual int get_ex_overlap () {
			return 1 + n_boundary_out;
		}
		
		virtual int get_excess () {
			return n_boundary_in;
		}
		
		virtual void send (datatype *data_temp, int lda) {
			buffer.resize ((std::max (n_boundary_out, n_boundary_in) + 1) * ldn);
			for (int i = 0; i < ldn; ++i) {
				for (int j = 0; j < n_boundary_out + 1; ++j) {
					buffer [j * ldn + i] = data_temp [i * lda + j];
				}
			}

			messenger_ptr->send ((n_boundary_out + 1) * ldn, &buffer [0], out_id, 0);
		}
		
		virtual void receive (datatype *data_temp, int lda) {
			buffer.resize ((std::max (n_boundary_out, n_boundary_in) + 1) * ldn);
			messenger_ptr->recv ((n_boundary_in + 1) * ldn, &buffer [0], out_id, 0);

			for (int i = 0; i < ldn; ++i) {
				for (int j = 0; j < n_boundary_out + 1; ++j) {
					data_temp [i * lda + j] += buffer [j * ldn + i];
				}
			}
		}
		
		virtual void calculate_rhs (datatype *data, datatype *interpolate_original, datatype *interpolate_data, datatype *data_temp, int m, int lda, int flag) {
			linalg::matrix_scale (1 + n_boundary_out + n_boundary_in, ldn, 0.0, data_temp + (top ? 1 : -1 - n_boundary_out - n_boundary_in), lda);
			// Setting the external overlapping boundary
			// linalg::interpolate (n_boundary_out, ldn, m, 1.0, 1.0, positions, interpolate_data, boundary_positions, data_temp + (top ? 1 + n_boundary_in : -n_boundary_in - n_boundary_out), lda, lda);
			// Scale the internal boundary row
			linalg::scale (ldn, alpha, data_temp, lda);
			// Copy the internal boundary row to the external boundary row
			linalg::copy (ldn, data_temp, data_temp + (top ? 1 : -1) * (1 + n_boundary_out + n_boundary_in), lda, lda);
			// Add the original data to the overlapping boundary
			// if (interpolate_original) {
			// 	linalg::interpolate (n_boundary_out, ldn, m, 1.0, 1.0, positions, interpolate_original, boundary_positions, data_temp + (top ? 1 + n_boundary_in : -n_boundary_in - n_boundary_out), m, lda);
			// }
			// Add the original data to the internal boundary row
			if (data) {
				linalg::add_scaled (ldn, 1.0, data, data_temp, m, lda);
			}
		}
		
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *interpolate_matrix, datatype *matrix_out, int lda, bool diverging = false) {
			// Zero everything but the internal boundary row
			linalg::matrix_scale (1 + n_boundary_out + n_boundary_in, m, 0.0, matrix_out + (top ? 1 : -1 - n_boundary_out - n_boundary_in), lda);
			// Setting the external boundary matrix row
			linalg::matrix_add_scaled (1, m, alpha, matrix_in, matrix_out + (top ? 1 : -1) * (1 + n_boundary_out + n_boundary_in), m, lda);
			// Setting the external overlapping boundary matrix row
			// linalg::interpolate (n_boundary_out, m, m, 1.0, 1.0, positions, interpolate_matrix, boundary_positions, matrix_out + (top ? 1 + n_boundary_in : -n_boundary_in - n_boundary_out), m, lda);
			// Setting the internal boundary matrix row
			linalg::matrix_add_scaled (1, m, alpha, matrix_in, matrix_out, m, lda);
			linalg::matrix_scale (1 + n_boundary_out + n_boundary_in + 1, m, timestep, matrix_out + (top ? 0 : -2 - n_boundary_out - n_boundary_in), lda);
			linalg::matrix_add_scaled (1 + n_boundary_in, m, 1.0, default_matrix + (top ? 0 : -n_boundary_in), matrix_out + (top ? 0 : -n_boundary_in), m, lda);
			if (diverging) {
				linalg::matrix_add_scaled (1, m, -1.0, default_matrix, matrix_out + (top ? 1 + n_boundary_in + n_boundary_out : -n_boundary_in - n_boundary_out - 1), m, lda);
			}
			linalg::interpolate (n_boundary_out, m, m, -1.0, 0.0, positions, interpolate_matrix, boundary_positions, matrix_out + (top ? 1 + n_boundary_in : -n_boundary_in - n_boundary_out), m, lda);
			
			DEBUG ("Calculated.");
		}
	};
} /* plans */

#endif /* end of include guard: BOUNDARY_TWO_D_HPP_863EDB36 */