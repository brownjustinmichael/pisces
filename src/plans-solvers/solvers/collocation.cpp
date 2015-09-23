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
	namespace solvers
	{
		template <class datatype>
		void collocation <datatype>::init (mpi::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype *i_rhs, grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out) {
			TRACE ("Building solver...");
			n = i_data.get_grid (0).get_n ();
			ldn = i_data.get_grid (0).get_ld ();
			m = i_data.get_grid (1).get_n ();
			messenger_ptr = i_messenger_ptr;
			excess_0 = i_data.get_grid (1).get_excess_0 ();
			excess_n = i_data.get_grid (1).get_excess_n ();
			default_matrix = i_data.get_grid (1).get_data (0);

			matrix.resize (m * m, 0.0);
			rhs_ptr = i_rhs;
			
			boundary_0 = i_boundary_0;
			boundary_n = i_boundary_n;
			
			// For convenience and speed, get the overlap information now
			ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
			overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
			ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
			overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
			lda = m + ex_overlap_n + ex_overlap_0;
			inner_m = lda - overlap_0 - overlap_n;
			int ns0 = overlap_0;
			
			// Gather the extents of all the elements, and resize the boundary matrix information to fit
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
			
			// Resize the matrices and right hand sides for the data
			factorized_matrix.resize (lda * lda, 0.0);
			ipiv.resize (m); // TODO Should be n - ntop - nbot - excess_0 - excess_n
			data_temp.resize (lda * ldn);
			TRACE ("Solver built.");
		}
	
		template <class datatype>
		void collocation <datatype>::factorize () {
			int info;
			std::stringstream debug;
			
			TRACE ("Factorizing...");
			
			// Copy over the left hand side into the matrix
			linalg::matrix_copy (m, m, &matrix [0], &factorized_matrix [(ex_overlap_0) * (lda + 1)], m, lda);
			
			/*
				TODO Should we do the matrix copy before the edges?
			*/
			
			// Calculate the matrix contributions from the boundaries. Note, these take care of the timestep multiplication
			if (boundary_0) {
				boundary_0->calculate_matrix (timestep / 2.0, default_matrix + excess_0, &matrix [excess_0], &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0], lda);
			}
			if (boundary_n) {
				boundary_n->calculate_matrix (timestep / 2.0, default_matrix + m - 1 - excess_n, &matrix [m - 1 - excess_n], &factorized_matrix [(ex_overlap_0) * (lda + 1) + m - 1 - excess_n], lda);
			}
			
			// Scale the matrices by the timestep
			linalg::matrix_scale (lda - 2 - excess_0 - ex_overlap_0 - excess_n - ex_overlap_n, lda, timestep / 2.0, &factorized_matrix [ex_overlap_0 + 1 + excess_0], lda);
			
			// Add in the original data to the matrix
			linalg::matrix_add_scaled (m - 2 - excess_0 - excess_n, m, 1.0, default_matrix + excess_0 + 1, &factorized_matrix [(ex_overlap_0) * (lda + 1) + excess_0 + 1], m, lda);
			
			linalg::block::matrix_factorize (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, lda, sqrt ((int) boundary_matrix.size ()));
			
			TRACE ("Done.");
		}
	
		template <class datatype>
		void collocation <datatype>::execute () {
			solver <datatype>::execute ();
			int info;
			TRACE ("Executing solve...");
			
			/*
				TODO Add timestep check here?
			*/
			
			linalg::scale ((ldn) * lda, 0.0, &data_temp [0]);
			
			// Add the right hand side, scaled by the timestep
			if (rhs_ptr) {
				linalg::matrix_add_scaled (m, ldn, timestep / 2.0, rhs_ptr, &data_temp [ex_overlap_0], m, lda);
			}
			
			// Include the contributions from the boundaries; if there aren't boundaries, just use the data present
			if (boundary_0) {
				boundary_0->calculate_rhs (data_in + excess_0, &data_temp [ex_overlap_0 + excess_0], m, lda, z_solve);
			} else {
				linalg::matrix_add_scaled (1 + excess_0, ldn, 1.0, data_in, &data_temp [ex_overlap_0], m, lda);
			}
			if (boundary_n) {
				boundary_n->calculate_rhs (data_in + m - 1 - excess_n, &data_temp [lda - 1 - excess_n - ex_overlap_n], m, lda, z_solve);
			} else {
				linalg::matrix_add_scaled (1 + excess_n, ldn, 1.0, data_in + m - 1 - excess_n, &data_temp [lda - 1 - excess_n - ex_overlap_n], m, lda);
			}
			
			// Add the data from the previous timestep
			linalg::matrix_add_scaled (m - 2 - excess_0 - excess_n, ldn, 1.0, data_in + 1 + excess_0, &data_temp [ex_overlap_0 + 1 + excess_0], m, lda);
			
			TRACE ("Performing matrix solve.")

			linalg::block::matrix_solve (messenger_ptr->get_id (), messenger_ptr->get_np (), inner_m, overlap_0, overlap_n, &factorized_matrix [0], &ipiv [0], &data_temp [0], &boundary_matrix [0], messenger_ptr->get_id () == 0 ? &bipiv [0] : NULL, messenger_ptr->get_id () == 0 ? &ns [0] : NULL, &info, ldn, lda, sqrt ((int) boundary_matrix.size ()), lda);
			
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
			linalg::matrix_copy (m, ldn, &data_temp [ex_overlap_0], data_out, lda, m);
			
			// This solve transforms, so we specify that in the component flags
			component_flags |= transformed_vertical;
			
			TRACE ("Solve complete.")
		}
		
		template class collocation <double>;
	} /* solvers */
} /* plans */

#endif /* end of include guard: COLLOCATION_CPP_655F168D */
