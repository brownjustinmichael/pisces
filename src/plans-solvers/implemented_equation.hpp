/*!**********************************************************************
 * \file implemented_equation.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IMPLEMENTED_EQUATION_HPP_74410EFA
#define IMPLEMENTED_EQUATION_HPP_74410EFA

#include "equation.hpp"
#include "plans-transforms/transform.hpp"

namespace plans
{
	template <class datatype>
	class implemented_equation : public plans::equation <datatype>
	{
	private:
		int n;
		int ldn;
		int m;
		int flags;
		plans::grid <datatype> &grid_n;
		plans::grid <datatype> &grid_m;
		
		std::shared_ptr <plans::solver <datatype> > x_solver;
		std::shared_ptr <plans::solver <datatype> > z_solver;
		
		std::vector <datatype> spectral_rhs_vec;
		std::vector <datatype> real_rhs_vec;
		std::vector <datatype> old_rhs_vec;
		std::vector <datatype> new_rhs_vec;
		std::vector <datatype> cor_rhs_vec;
		
		datatype *spectral_rhs_ptr;
		datatype *real_rhs_ptr, *old_rhs_ptr, *new_rhs_ptr, *cor_rhs_ptr;
		
		std::shared_ptr <plans::plan <datatype> > transform;
		
		using plans::equation <datatype>::data;
		using plans::equation <datatype>::element_flags;
		using plans::equation <datatype>::component_flags;
		
	public:
		implemented_equation (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype *i_data, int *i_element_flags, int *i_component_flags) : plans::equation <datatype> (i_data, i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), grid_n (i_grid_n), grid_m (i_grid_m) {
			spectral_rhs_ptr = NULL;
			real_rhs_ptr = NULL;
			old_rhs_vec.resize (ldn * m, 0.0);
			old_rhs_ptr = &old_rhs_vec [0];
			new_rhs_vec.resize (ldn * m, 0.0);
			new_rhs_ptr = &new_rhs_vec [0];
			cor_rhs_vec.resize (ldn * m, 0.0);
			cor_rhs_ptr = &cor_rhs_vec [0];
		}
		
		virtual ~implemented_equation () {}
		
		virtual int n_dependencies () {
			if (*component_flags & x_solve) {
				if (x_solver) {
					return x_solver->n_dependencies ();
				} else {
					return 0;
				}
			} else if (*component_flags & z_solve) {
				if (z_solver) {
					return z_solver->n_dependencies ();
				} else {
					return 0;
				}
			} else {
				return 0;
			}
		}

		virtual const std::string& get_dependency (int i) {
			if (*component_flags & x_solve) {
				return x_solver->get_dependency (i);
			} else if (*component_flags & z_solve) {
				return z_solver->get_dependency (i);
			} else {
				FATAL ("Cannot get dependency without valid solve direction " << *component_flags);
				throw 0;
			}
		}

		virtual void add_dependency (std::string name, int flags = 0x00) {
			if (!(flags & not_x_solver)) {
				if (x_solver) {
					x_solver->add_dependency (name);
				}
			} else if (!(flags & not_z_solver)) {
				if (z_solver) {
					z_solver->add_dependency (name);
				}
			}
		}
		
		plans::grid <datatype> *grid_ptr (int index = 0) {
			if (index == 0) {
				return &grid_n;
			} else {
				return &grid_m;
			}
		}
		
		datatype *rhs_ptr (int index = spectral_rhs) {
			if (index == spectral_rhs) {
				if (!spectral_rhs_ptr) {
					spectral_rhs_vec.resize (ldn * m);
					spectral_rhs_ptr = &spectral_rhs_vec [0];
				}
				return spectral_rhs_ptr;
			} else if (index == real_rhs) {
				if (!real_rhs_ptr) {
					real_rhs_vec.resize (ldn * m);
					real_rhs_ptr = &real_rhs_vec [0];
					flags = 0x00;
					transform = std::shared_ptr <plans::plan <datatype> > (new plans::horizontal_transform <datatype> (n, m, real_rhs_ptr, NULL, 0x00, element_flags, &flags));
				}
				return real_rhs_ptr;
			} else {
				return NULL;
			}
		}
		
		datatype *matrix_ptr (int index = 0) {
			if (index == 0) {
				if (x_solver) {
					return x_solver->matrix_ptr ();
				} else {
					return NULL;
				}
			} else {
				if (z_solver) {
					return z_solver->matrix_ptr ();
				} else {
					return NULL;
				}
			}
		}
		
		virtual void reset () {
			if (spectral_rhs_ptr) {
				linalg::scale (ldn * m, 0.0, spectral_rhs_ptr);
			}
			if (real_rhs_ptr) {
				linalg::scale (ldn * m, 0.0, real_rhs_ptr);
			}
			linalg::scale (m * ldn, 0.0, new_rhs_ptr);
			
			if (*component_flags & z_solve) {
				*component_flags &= ~z_solve;
				*component_flags |= x_solve;
			} else {
				*component_flags &= ~x_solve;
				*component_flags |= z_solve;
			}
		}
		
		virtual void add_solver (std::shared_ptr <plans::solver <datatype>> i_solver, int flags = 0x00) {
			TRACE ("Adding solver...");
			if (!(flags & not_x_solver)) {
				x_solver = i_solver;
			}
			if (!(flags & not_z_solver)) {
				z_solver = i_solver;
			}
			TRACE ("Added.");
		}
		
		virtual void add_solver (const typename plans::solver <datatype>::factory &factory, int flags = 0x00) {
			plans::grid <datatype>* grids [2] = {&grid_n, &grid_m};
			plans::implemented_equation <datatype>::add_solver (factory.instance (grids, data, cor_rhs_ptr, element_flags, component_flags), flags);
		}
		
		virtual std::shared_ptr <plans::solver <datatype>> get_solver (int flags = 0x00) {
			if (!(flags & not_x_solver)) {
				return x_solver;
			}
			if (!(flags & not_z_solver)) {
				return z_solver;
			}
			FATAL ("Invalid flags: " << flags);
			throw 0;
		}
		
		void add_plan (const typename plans::explicit_plan <datatype>::factory &factory, int flags) {
			TRACE ("Adding plan...");
			plans::grid <datatype>* grids [2] = {&grid_n, &grid_m};
			plans::equation <datatype>::add_plan (factory.instance (grids, data, new_rhs_ptr, element_flags, component_flags), flags);
		}
		
		void add_plan (const typename plans::real_plan <datatype>::factory &factory, int flags) {
			TRACE ("Adding plan...");
			plans::grid <datatype>* grids [2] = {&grid_n, &grid_m};
			plans::equation <datatype>::add_plan (factory.instance (grids, data, rhs_ptr (real_rhs), element_flags, component_flags), flags);
		}
		
		void add_plan (const typename plans::implicit_plan <datatype>::factory &factory, int flags) {
			TRACE ("Adding plan...");
			plans::grid <datatype>* grids [2] = {&grid_n, &grid_m};
			datatype* matrices [2] = {matrix_ptr (0), matrix_ptr (1)};
			plans::equation <datatype>::add_plan (factory.instance (grids, matrices, data, rhs_ptr (spectral_rhs), element_flags, component_flags), flags);
		}
		
	protected:
		virtual void _factorize () {
			if (x_solver) {
				x_solver->factorize ();
			}
			if (z_solver && (x_solver != z_solver)) {
				z_solver->factorize ();
			}
		}
		
		virtual void _solve () {
			TRACE ("Solving...");
			if (transform) transform->execute ();
			
			// linalg::matrix_scale (m, ldn, 0.0, cor_rhs_ptr);
			//
			// // De-alias the RHS
			// if (real_rhs_ptr) {
			// 	linalg::matrix_scale (m, ldn / 3, 0.0, real_rhs_ptr + m * (ldn - ldn / 3));
			// 	linalg::matrix_copy (m, ldn, real_rhs_ptr, cor_rhs_ptr);
			// }
			//
			// if (spectral_rhs_ptr) {
			// 	linalg::matrix_add_scaled (m, ldn, 1.0, real_rhs_ptr, cor_rhs_ptr);
			// 	linalg::matrix_scale (m, ldn, -0.0, old_rhs_ptr);
			// 	linalg::matrix_add_scaled (m, ldn, 1.0, cor_rhs_ptr, cor_rhs_ptr);
			// 	// linalg::matrix_copy (m, ldn, spectral_rhs_ptr, old_rhs_ptr);
			// }
			
			linalg::matrix_copy (m, ldn, old_rhs_ptr, cor_rhs_ptr);
			linalg::matrix_scale (m, ldn, -0.5, cor_rhs_ptr); // Set to -0.5 for AB

			// De-alias the RHS
			if (real_rhs_ptr) {
				linalg::matrix_scale (m, ldn / 3, 0.0, real_rhs_ptr + m * (ldn - ldn / 3));
				linalg::matrix_add_scaled (m, ldn, 1.0, real_rhs_ptr, new_rhs_ptr);
			}
		
			linalg::matrix_add_scaled (m, ldn, 1.5, new_rhs_ptr, cor_rhs_ptr); // Set to 1.5 for AB
			linalg::matrix_copy (m, ldn, new_rhs_ptr, old_rhs_ptr);

			if (spectral_rhs_ptr) {
				linalg::matrix_add_scaled (m, ldn, 1.0, spectral_rhs_ptr, cor_rhs_ptr);
			}
			
			if (*component_flags & x_solve) {
				if (x_solver) {
					x_solver->execute ();
				}
			} else if (*component_flags & z_solve) {
				if (z_solver) {
					z_solver->execute ();
				}
			}
		}
	};
} /* plans */

#endif /* end of include guard: IMPLEMENTED_EQUATION_HPP_74410EFA */
