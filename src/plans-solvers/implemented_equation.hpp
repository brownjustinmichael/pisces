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
	namespace solvers
	{
		/*!**********************************************************************
		 * \brief An implementation of the equation class in 2D
		 * 
		 * \copydoc equation
		 * 
		 * This implements the equation class into a two-dimensional object. It should be possible to remove this class by making equation more general.
		 ************************************************************************/
		template <class datatype>
		class implemented_equation : public equation <datatype>
		{
		private:
			using plans::solvers::equation <datatype>::data;
			using plans::solvers::equation <datatype>::element_flags;
			using plans::solvers::equation <datatype>::component_flags;
			
			int n; //!< The integer number of elements in the horizontal
			int ldn; //!< The integer length of the array in the horizontal (can be larger than n)
			int m; //!< The integer number of elements in the vertical
			int flags; //!< The integer flags associated with the equation
			grids::grid <datatype> &grid_n; //!< The grid object in the horizontal
			grids::grid <datatype> &grid_m; //!< The grid object in the vertical
			
			std::shared_ptr <plans::solvers::solver <datatype> > x_solver; //!< A pointer to the horizontal solver
			std::shared_ptr <plans::solvers::solver <datatype> > z_solver; //!< A pointer to the vertical solver
			
			std::vector <datatype> spectral_rhs_vec; //!< A vector containing the right hand side from the spectral terms
			std::vector <datatype> real_rhs_vec; //!< A vector containing the right hand side from the real terms
			std::vector <datatype> old_rhs_vec; //!< A vector containing the right hand side from the last timestep
			std::vector <datatype> old2_rhs_vec; //!< A vector containing the right hand side from two timesteps ago
			std::vector <datatype> old3_rhs_vec; //!< A vector containing the right hand side from three timesteps ago
			std::vector <datatype> new_rhs_vec; //!< A vector containing the total current right hand side
			std::vector <datatype> cor_rhs_vec; //!< A vector containing the right hand side after the Adams Bashforth scheme
			
			datatype *spectral_rhs_ptr; //!< A pointer to the spectral_rhs_vec for speed
			datatype *real_rhs_ptr; //!< A pointer to the real_rhs_vec for speed
			datatype *old_rhs_ptr; //!< A pointer to the old_rhs_vec for speed
			datatype *old2_rhs_ptr; //!< A pointer to the old2_rhs_vec for speed
			datatype *old3_rhs_ptr; //!< A pointer to the old3_rhs_vec for speed
			datatype *new_rhs_ptr; //!< A pointer to the new_rhs_vec for speed
			datatype *cor_rhs_ptr; //!< A pointer to the cor_rhs_vec for speed
			
			std::shared_ptr <plans::plan <datatype> > transform; //!< A shared pointer to a transform if the equation has a real_rhs_vec
		
		public:
			/*!**********************************************************************
			 * \copydoc equation::equation
			 * 
			 * \param i_grid_n The grid object in the horizontal
			 * \param i_grid_m The grid object in the vertical
			 ************************************************************************/
			implemented_equation (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype *i_data, int *i_element_flags, int *i_component_flags) : plans::solvers::equation <datatype> (i_data, i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), grid_n (i_grid_n), grid_m (i_grid_m) {
				spectral_rhs_ptr = NULL;
				real_rhs_ptr = NULL;
				old_rhs_vec.resize (ldn * m, 0.0);
				old_rhs_ptr = &old_rhs_vec [0];
				old2_rhs_vec.resize (ldn * m, 0.0);
				old2_rhs_ptr = &old2_rhs_vec [0];
				old3_rhs_vec.resize (ldn * m, 0.0);
				old3_rhs_ptr = &old3_rhs_vec [0];
				new_rhs_vec.resize (ldn * m, 0.0);
				new_rhs_ptr = &new_rhs_vec [0];
				cor_rhs_vec.resize (ldn * m, 0.0);
				cor_rhs_ptr = &cor_rhs_vec [0];
			}
			
			virtual ~implemented_equation () {}
			
			/*!**********************************************************************
			 * \copydoc equation::n_dependencies
			 ************************************************************************/
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
			
			/*!**********************************************************************
			 * \copydoc equation::get_dependency
			 ************************************************************************/
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
			
			/*!**********************************************************************
			 * \copydoc equation::add_dependency
			 ************************************************************************/
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
			
			/*!**********************************************************************
			 * \copydoc equation::grid_ptr
			 ************************************************************************/
			grids::grid <datatype> *grid_ptr (int index = 0) {
				if (index == 0) {
					return &grid_n;
				} else {
					return &grid_m;
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::rhs_ptr
			 ************************************************************************/
			datatype *rhs_ptr (int index = spectral_rhs) {
				if (index == spectral_rhs) {
					// If the spectral_rhs has been requested and doesn't exist, make it
					if (!spectral_rhs_ptr) {
						spectral_rhs_vec.resize (ldn * m);
						spectral_rhs_ptr = &spectral_rhs_vec [0];
					}
					return spectral_rhs_ptr;
				} else if (index == real_rhs) {
					// If the real_rhs has been requested and doesn't exist, make it and its transform
					if (!real_rhs_ptr) {
						real_rhs_vec.resize (ldn * m);
						real_rhs_ptr = &real_rhs_vec [0];
						flags = 0x00;
						transform = std::shared_ptr <plans::plan <datatype> > (new plans::transforms::horizontal <datatype> (n, m, real_rhs_ptr, NULL, 0x00, element_flags, &flags));
					}
					return real_rhs_ptr;
				} else {
					// Use the default rhs
					return new_rhs_ptr;
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::matrix_ptr
			 ************************************************************************/
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
			
			/*!**********************************************************************
			 * \copydoc equation::reset
			 ************************************************************************/
			virtual void reset () {
				if (spectral_rhs_ptr) {
					linalg::scale (ldn * m, 0.0, spectral_rhs_ptr);
				}
				if (real_rhs_ptr) {
					linalg::scale (ldn * m, 0.0, real_rhs_ptr);
				}
				linalg::scale (m * ldn, 0.0, new_rhs_ptr);
				
				// Since this is an alternating direction solve, make sure to switch directions on reset
				if (*component_flags & z_solve) {
					*component_flags &= ~z_solve;
					*component_flags |= x_solve;
				} else {
					*component_flags &= ~x_solve;
					*component_flags |= z_solve;
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::get_solver
			 ************************************************************************/
			virtual std::shared_ptr <plans::solvers::solver <datatype>> get_solver (int flags = 0x00) {
				if (!(flags & not_x_solver)) {
					return x_solver;
				}
				if (!(flags & not_z_solver)) {
					return z_solver;
				}
				FATAL ("Invalid flags: " << flags);
				throw 0;
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_solver
			 ************************************************************************/
			virtual void add_solver (std::shared_ptr <plans::solvers::solver <datatype>> i_solver, int flags = 0x00) {
				TRACE ("Adding solver...");
				if (!(flags & not_x_solver)) {
					x_solver = i_solver;
				}
				if (!(flags & not_z_solver)) {
					z_solver = i_solver;
				}
				TRACE ("Added.");
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_solver(const typename solver<datatype>::factory&,int)
			 ************************************************************************/
			virtual void add_solver (const typename plans::solvers::solver <datatype>::factory &i_factory, int flags = 0x00) {
				grids::grid <datatype>* grids [2] = {&grid_n, &grid_m};
				plans::solvers::implemented_equation <datatype>::add_solver (i_factory.instance (grids, data, cor_rhs_ptr, element_flags, component_flags), flags);
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_plan(const typename explicit_plan<datatype>::factory&,int)
			 ************************************************************************/
			void add_plan (const typename plans::explicit_plan <datatype>::factory &i_factory, int flags) {
				TRACE ("Adding plan...");
				grids::grid <datatype>* grids [2] = {&grid_n, &grid_m};
				plans::solvers::equation <datatype>::add_plan (i_factory.instance (grids, data, new_rhs_ptr, element_flags, component_flags), flags);
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_plan(const typename explicit_plan<datatype>::factory&,int)
			 ************************************************************************/
			void add_plan (const typename plans::real_plan <datatype>::factory &i_factory, int flags) {
				TRACE ("Adding plan...");
				grids::grid <datatype>* grids [2] = {&grid_n, &grid_m};
				plans::solvers::equation <datatype>::add_plan (i_factory.instance (grids, data, rhs_ptr (real_rhs), element_flags, component_flags), flags);
			}
			
			/*!**********************************************************************
			 * \copydoc equation::add_plan(const typename explicit_plan<datatype>::factory&,int)
			 ************************************************************************/
			void add_plan (const typename plans::implicit_plan <datatype>::factory &i_factory, int flags) {
				TRACE ("Adding plan...");
				grids::grid <datatype>* grids [2] = {&grid_n, &grid_m};
				datatype* matrices [2] = {matrix_ptr (0), matrix_ptr (1)};
				plans::solvers::equation <datatype>::add_plan (i_factory.instance (grids, matrices, data, rhs_ptr (spectral_rhs), element_flags, component_flags), flags);
			}
		
		protected:
			/*!**********************************************************************
			 * \copydoc equation::_factorize
			 ************************************************************************/
			virtual void _factorize () {
				if (x_solver) {
					x_solver->factorize ();
				}
				if (z_solver && (x_solver != z_solver)) {
					z_solver->factorize ();
				}
			}
			
			/*!**********************************************************************
			 * \copydoc equation::_solve
			 ************************************************************************/
			virtual void _solve () {
				TRACE ("Solving...");
				
				// Transform the real_rhs
				if (transform) transform->execute ();
				
				// Add in the components from the last three timesteps for the AB scheme
				linalg::matrix_copy (m, ldn, old3_rhs_ptr, cor_rhs_ptr);
				linalg::matrix_scale (m, ldn, -3. / 8., cor_rhs_ptr);
				linalg::matrix_add_scaled (m, ldn, 37. / 24., old2_rhs_ptr, cor_rhs_ptr);
				linalg::matrix_add_scaled (m, ldn, -59. / 24., old_rhs_ptr, cor_rhs_ptr);
				
				// De-alias the RHS
				if (real_rhs_ptr) {
					linalg::matrix_scale (m, ldn / 3, 0.0, real_rhs_ptr + m * (ldn - ldn / 3));
					linalg::matrix_add_scaled (m, ldn, 1.0, real_rhs_ptr, new_rhs_ptr);
				}
				
				// Add in the component from this timestep
				linalg::matrix_add_scaled (m, ldn, 55. / 24., new_rhs_ptr, cor_rhs_ptr);
				
				// Rotate the old timestep pointers
				datatype *temp = old3_rhs_ptr;
				old3_rhs_ptr = old2_rhs_ptr;
				old2_rhs_ptr = old_rhs_ptr;
				old_rhs_ptr = temp;
				linalg::matrix_copy (m, ldn, new_rhs_ptr, old_rhs_ptr);
				
				// Add in the spectral component (the explicit part of the implicit component) after the AB scheme
				if (spectral_rhs_ptr) {
					linalg::matrix_add_scaled (m, ldn, 1.0, spectral_rhs_ptr, cor_rhs_ptr);
				}
				
				// Solve either the x direction solve or the z direction solve
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
	} /* solvers */
} /* plans */

#endif /* end of include guard: IMPLEMENTED_EQUATION_HPP_74410EFA */
