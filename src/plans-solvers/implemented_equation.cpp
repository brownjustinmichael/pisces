/*!**********************************************************************
 * \file implemented_equation.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include "implemented_equation.hpp"

namespace plans
{
	namespace solvers
	{
		void implemented_equation::_solve () {
			TRACE ("Solving...");
			static int count = 0;
			
			// Transform the real_rhs
			if (transform) transform->execute ();
			
			// Add in the components from the last three timesteps for the AB scheme
			linalg::matrix_copy (m, ldn, old3_rhs_ptr->ptr (real_spectral), cor_rhs_ptr->ptr (real_spectral));
			linalg::matrix_scale (m, ldn, -3. / 8., cor_rhs_ptr->ptr (real_spectral));
			linalg::matrix_add_scaled (m, ldn, 37. / 24., old2_rhs_ptr->ptr (real_spectral), cor_rhs_ptr->ptr (real_spectral));
			linalg::matrix_add_scaled (m, ldn, -59. / 24., old_rhs_ptr->ptr (real_spectral), cor_rhs_ptr->ptr (real_spectral));
			
			// De-alias the RHS
			linalg::matrix_scale (m, ldn / 3, 0.0, new_rhs_ptr->ptr (real_real) + m * (ldn - ldn / 3));
			linalg::matrix_add_scaled (m, ldn, 1.0, new_rhs_ptr->ptr (real_real), new_rhs_ptr->ptr (real_spectral));
			
			// Add in the component from this timestep
			linalg::matrix_add_scaled (m, ldn, 55. / 24., new_rhs_ptr->ptr (real_spectral), cor_rhs_ptr->ptr (real_spectral));

			// // Testing first order
			// linalg::matrix_copy (m, ldn, new_rhs_ptr->ptr (real_spectral), cor_rhs_ptr->ptr (real_spectral));
			
			// Rotate the old timestep pointers
			std::shared_ptr <grids::variable> temp = old3_rhs_ptr;
			old3_rhs_ptr = old2_rhs_ptr;
			old2_rhs_ptr = old_rhs_ptr;
			old_rhs_ptr = temp;
			linalg::matrix_copy (m, ldn, new_rhs_ptr->ptr (real_spectral), old_rhs_ptr->ptr (real_spectral));
			
			// Add in the spectral component (the explicit part of the implicit component) after the AB scheme
			// linalg::matrix_add_scaled (m, ldn, 1.0, new_rhs_ptr->ptr (spectral_spectral), cor_rhs_ptr->ptr (real_spectral));
			if (*component_flags & ignore_net) {
				linalg::scale (2 * m, 0.0, cor_rhs_ptr->ptr (real_spectral));
			}
			
			int state = -1;
			// Solve either the x direction solve or the z direction solve
			if ((*component_flags & x_solve)) {
				if (x_solver) {
					x_solver->execute ();
					state = x_solver->get_state ();
				}

				if (z_solver) {
					if (state >= 0 && state != z_solver->get_state_in ()) {
						transformer->update ();
					}

					linalg::matrix_scale (m, ldn, 0.0, cor_rhs_ptr->ptr (real_spectral));

					z_solver->execute ();
				}

				*component_flags &= ~x_solve;
				*component_flags |= z_solve;
			} else {
				// Solve either the x direction solve or the z direction solve
				if (z_solver) {
					z_solver->execute ();
					state = z_solver->get_state ();
				}

				if (x_solver) {
					if (state >= 0 && state != x_solver->get_state_in ()) {
						transformer->update ();
					}

					linalg::matrix_scale (m, ldn, 0.0, cor_rhs_ptr->ptr (real_spectral));

					x_solver->execute ();

					*component_flags &= ~z_solve;
					*component_flags |= x_solve;
				}
			}
		}
	} /* solvers */
} /* plans */
