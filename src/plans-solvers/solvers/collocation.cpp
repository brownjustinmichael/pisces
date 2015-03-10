/*!**********************************************************************
 * \file collocation.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_CPP_655F168D
#define COLLOCATION_CPP_655F168D

#include "collocation.hpp"

#include "linalg/exceptions.hpp"

namespace plans
{
	template <class datatype>
	collocation_solver <datatype>::collocation_solver (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
	plans::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), messenger_ptr (i_messenger_ptr), timestep (i_timestep), positions (&(i_grid_m [0])), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()), default_matrix (i_grid_m.get_data (0)) {
		TRACE ("Building solver...");
		matrix.resize (m * m, 0.0);
		rhs_ptr = i_rhs;
		
		boundary_0 = i_boundary_0;
		boundary_n = i_boundary_n;
		
		ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
		overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
		ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
		overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
		lda = m + ex_overlap_n + ex_overlap_0;
		inner_m = lda - overlap_0 - overlap_n;
		int ns0 = overlap_0;
		
		if (messenger_ptr->get_id () == 0) {
			ns.resize (messenger_ptr->get_np ());
			messenger_ptr->template gather <int> (1, &ns0, &ns [0]);
			int ntot = 0;
			for (int i = 0; i < messenger_ptr->get_np (); ++i) {
				ntot += ns [i];
			}
			boundary_matrix.resize (ntot * ntot);
			bipiv.resize (ntot);
		} else {
			messenger_ptr->template gather <int> (1, &ns0, NULL);
			boundary_matrix.resize ((overlap_0 + overlap_n) * (overlap_0 + overlap_n));
		}
		
		factorized_matrix.resize (lda * lda, 0.0);
		ipiv.resize (m); // TODO Should be n - ntop - nbot - excess_0 - excess_n
		data_temp.resize (lda * ldn);
		TRACE ("Solver built.");
	}
	
	template <class datatype>
	void collocation_solver <datatype>::factorize () {
		int info;
		std::stringstream debug;

		TRACE ("Factorizing...");

		// linalg::matrix_scale (lda, lda, 0.0, &factorized_matrix [0]);
		linalg::matrix_copy (m, m, &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);

		/*
			TODO Should we do the matrix copy before the edges?
		*/

		// DEBUG ("ZERO POINT " << &factorized_matrix [0]);

		if (boundary_0) {
			boundary_0->calculate_matrix (timestep, default_matrix + excess_0, &matrix [excess_0], &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0], lda);
		}
		if (boundary_n) {
			boundary_n->calculate_matrix (timestep, default_matrix + m - 1 - excess_n, &matrix [m - 1 - excess_n], &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1) + m - 1 - excess_n], lda);
		}
		
		linalg::matrix_scale (lda - 2 - excess_0 - ex_overlap_0 - excess_n - ex_overlap_n, lda, timestep, &factorized_matrix [ex_overlap_0 + 1 + excess_0], lda);
		
		linalg::matrix_add_scaled (m - 2 - excess_0 - excess_n, m, 1.0, default_matrix + excess_0 + 1, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0 + 1], m, lda);
		
		// for (int j = 0; j < lda; ++j) {
		// 	for (int i = 0; i < lda; ++i) {
		// 		debug << factorized_matrix [i * lda + j] << " ";
		// 	}
		// 	DEBUG (debug.str ());
		// 	debug.str ("");
		// }
		
		// assert (false);

		linalg::p_block_matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));

		TRACE ("Done.");
	}
	
	template <class datatype>
	void collocation_solver <datatype>::execute () {
		int info;
		std::stringstream debug;
		TRACE ("Executing solve...");

		/*
			TODO Add timestep check here?
		*/
		
		// for (int j = 0; j < m; ++j) {
		// 	for (int i = 0; i < ldn; ++i) {
		// 		debug << data [i * m + j] << " ";
		// 	}
		// 	DEBUG ("DATA " << debug.str ());
		// 	debug.str ("");
		// }

		linalg::scale ((ldn) * lda, 0.0, &data_temp [0]);

		if (rhs_ptr) {
			linalg::matrix_add_scaled (m, ldn, timestep, rhs_ptr, &data_temp [ex_overlap_0], m, lda);
		}
		
		// for (int j = 0; j < lda; ++j) {
		// 	for (int i = 0; i < ldn; ++i) {
		// 		debug << data_temp [i * lda + j] << " ";
		// 	}
		// 	DEBUG ("RHS " << debug.str ());
		// 	debug.str ("");
		// }

		if (boundary_0) {
			boundary_0->calculate_rhs (data + excess_0, data, &data_temp [ex_overlap_0], &data_temp [ex_overlap_0 + excess_0], m, lda, z_solve);
		} else {
			linalg::matrix_add_scaled (1 + excess_0, ldn, 1.0, data, &data_temp [ex_overlap_0], m, lda);
		}
		if (boundary_n) {
			boundary_n->calculate_rhs (data + m - 1 - excess_n, data, &data_temp [ex_overlap_0], &data_temp [lda - 1 - excess_n - ex_overlap_n], m, lda, z_solve);
		} else {
			linalg::matrix_add_scaled (1 + excess_n, ldn, 1.0, data + m - 1 - excess_n, &data_temp [lda - 1 - excess_n - ex_overlap_n], m, lda);
		}
		
		linalg::matrix_add_scaled (m - 2 - excess_0 - excess_n, ldn, 1.0, data + 1 + excess_0, &data_temp [ex_overlap_0 + 1 + excess_0], m, lda);
	
		// for (int j = 0; j < lda; ++j) {
		// 	for (int i = 0; i < ldn; ++i) {
		// 		debug << data_temp [i * lda + j] << " ";
		// 	}
		// 	DEBUG ("RHS " << debug.str ());
		// 	debug.str ("");
		// }
		
		// assert (false);

		linalg::p_block_matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, ldn, lda, sqrt ((int) boundary_matrix.size ()), lda);

		// for (int j = 0; j < lda; ++j) {
		// 		for (int i = 0; i < ldn; ++i) {
		// 			debug << data_temp [i * lda + j] << " ";
		// 		}
		// 		DEBUG ("DONE " << debug.str ());
		// 		debug.str ("");
		// }

		TRACE ("Matrix solve complete.");

#ifdef CHECKNAN
		for (int i = 0; i < ldn; ++i) {
			for (int j = 0; j < m; ++j) {
				if (std::isnan (data_temp [ex_overlap_0 + i * m + j])) {
					FATAL ("Found nan.");
					for (int k = 0; k < m; ++k) {
						printf ("%f ", data_temp [ex_overlap_0 + k * m + j]);
					}
					printf ("\n");
					throw linalg::exceptions::nan ();
				}
			}
		}
#else
		if (std::isnan (data_temp [0])) {
			FATAL ("Found nan.");
			throw linalg::exceptions::nan ();
		}
#endif

		TRACE ("Updating...");
		linalg::matrix_copy (m, ldn, &data_temp [ex_overlap_0], data, lda, m);
		
		// for (int j = 0; j < m; ++j) {
		// 		for (int i = 0; i < ldn; ++i) {
		// 			debug << data [i * m + j] << " ";
		// 		}
		// 		DEBUG ("DONE " << debug.str ());
		// 		debug.str ("");
		// }
		
		*component_flags |= transformed_vertical;

		TRACE ("Solve complete.")
		TRACE ("Execution complete.");
	}
	
	template class collocation_solver <double>;
} /* plans */

#endif /* end of include guard: COLLOCATION_CPP_655F168D */
