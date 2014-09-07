/*!**********************************************************************
 * \file solver_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_TWO_D_HPP_JW4BV4PS
#define SOLVER_TWO_D_HPP_JW4BV4PS

#include "../utils/messenger.hpp"
#include "../bases/solver.hpp"
#include "transform_two_d.hpp"
#include "../utils/interpolate.hpp"
#include "plan_two_d.hpp"
#include "boundary_two_d.hpp"

namespace two_d
{
	template <class datatype>
	class master_solver : public bases::master_solver <datatype>
	{
	private:
		int n;
		int ldn;
		int m;
		int flags;
		bases::grid <datatype> &grid_n;
		bases::grid <datatype> &grid_m;
		
		std::shared_ptr <bases::solver <datatype> > x_solver;
		std::shared_ptr <bases::solver <datatype> > z_solver;
		
		std::vector <datatype> spectral_rhs_vec;
		std::vector <datatype> real_rhs_vec;
		
		datatype *spectral_rhs_ptr;
		datatype *real_rhs_ptr;
		
		std::shared_ptr <bases::plan <datatype> > transform;
		
		using bases::master_solver <datatype>::data;
		using bases::master_solver <datatype>::element_flags;
		using bases::master_solver <datatype>::component_flags;
		
	public:
		master_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_data, int *i_element_flags, int *i_component_flags) : bases::master_solver <datatype> (i_data, i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), grid_n (i_grid_n), grid_m (i_grid_m) {
			spectral_rhs_ptr = NULL;
			real_rhs_ptr = NULL;
		}
		
		virtual ~master_solver () {}
		
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

		virtual int& get_dependency (int i) {
			if (*component_flags & x_solve) {
				return x_solver->get_dependency (i);
			} else if (*component_flags & z_solve) {
				return z_solver->get_dependency (i);
			} else {
				throw 0;
			}
		}

		virtual void add_dependency (int name, int flags = 0x00) {
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
		
		bases::grid <datatype> *grid_ptr (int index = 0) {
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
					transform = std::shared_ptr <bases::plan <datatype> > (new fourier::horizontal_transform <datatype> (n, m, real_rhs_ptr, NULL, 0x00, element_flags, &flags));
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
				utils::scale (ldn * m, 0.0, spectral_rhs_ptr);
			}
			if (real_rhs_ptr) {
				utils::scale (ldn * m, 0.0, real_rhs_ptr);
			}
			
			if (*component_flags & z_solve) {
				*component_flags &= ~z_solve;
				*component_flags |= x_solve;
			} else {
				*component_flags &= ~x_solve;
				*component_flags |= z_solve;
			}
		}
		
		virtual void add_solver (std::shared_ptr <bases::solver <datatype>> i_solver, int flags = 0x00) {
			TRACE ("Adding solver...");
			if (!(flags & not_x_solver)) {
				DEBUG (1);
				x_solver = i_solver;
			}
			if (!(flags & not_z_solver)) {
				DEBUG (2);
				z_solver = i_solver;
			}
			TRACE ("Added.");
		}
		
		virtual std::shared_ptr <bases::solver <datatype>> get_solver (int flags = 0x00) {
			if (!(flags & not_x_solver)) {
				return x_solver;
			}
			if (!(flags & not_z_solver)) {
				return z_solver;
			}
			throw 0;
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
			
			if (spectral_rhs_ptr && real_rhs_ptr) {
				utils::matrix_add_scaled (m, ldn, 1.0, real_rhs_ptr, spectral_rhs_ptr);
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
	
	namespace fourier
	{
		template <class datatype>
		class collocation_solver : public bases::solver <datatype>
		{
		private:
			int n;
			int ldn;
			int m;
			datatype *data;
			int flags;

			utils::messenger* messenger_ptr;
	
			datatype& timestep; //!< A datatype reference to the current timestep
			datatype *rhs_ptr;

			const datatype* positions;
			int excess_0; //!< The integer number of elements to recv from edge_0
			int excess_n; //!< The integer number of elements to recv from edge_n

			datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component

			std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
			std::vector <datatype> boundary_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
			std::vector <datatype> previous_rhs;
			std::vector <int> ns;
			std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
			std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
			std::vector <datatype> matrix;
			
			std::shared_ptr <bases::boundary <datatype>> boundary_0, boundary_n;
			
			int inner_m;
			int ex_overlap_0;
			int overlap_0;
			int ex_overlap_n;
			int overlap_n;
			int lda;
			
			using bases::solver <datatype>::element_flags;
			using bases::solver <datatype>::component_flags;
			
		public:
			/*!**********************************************************************
			 * The collocation matrix is set up as 
			 * 
			 * 0 0 boundary row for above element       0 0
			 * 0 0 interpolating row for above element  0 0
			 * 0 0 [interpolating row for this element] 0 0
			 * 0 0 [boundary row for this element     ] 0 0
			 * 0 0 [matrix                            ] 0 0
			 * 0 0 [boundary row for this element     ] 0 0
			 * 0 0 [interpolating row for this element] 0 0
			 * 0 0 interpolating row for below element  0 0
			 * 0 0 boundary row for below element       0 0
			 ************************************************************************/
			collocation_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
			collocation_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n);
			
			virtual ~collocation_solver () {}
			
			datatype *matrix_ptr () {
				return &matrix [0];
			}
	
			void factorize ();
			void execute ();
		};
		
		template <class datatype>
		class fourier_solver : public bases::solver <datatype>
		{
		private:
			int n;
			int ldn;
			int m;
			datatype *data;
			int flags;
			
			datatype& timestep; //!< A datatype reference to the current timestep
			datatype *rhs_ptr;

			int excess_0; //!< The integer number of elements to recv from edge_0
			int excess_n; //!< The integer number of elements to recv from edge_n
			
			std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <datatype> previous_rhs;
			std::vector <datatype> horizontal_matrix;
			std::vector <datatype> factorized_horizontal_matrix;
			
			std::shared_ptr <bases::boundary <datatype>> boundary_0, boundary_n;
			
			int inner_m;
			int ex_overlap_0;
			int overlap_0;
			int ex_overlap_n;
			int overlap_n;
			int lda;
			const datatype *pos_m;
			
			using bases::solver <datatype>::element_flags;
			using bases::solver <datatype>::component_flags;
			
		public:
			/*!**********************************************************************
			 * The fourier matrix is set up as 
			 * 
			 * 0 0 boundary row for above element       0 0
			 * 0 0 interpolating row for above element  0 0
			 * 0 0 [interpolating row for this element] 0 0
			 * 0 0 [boundary row for this element     ] 0 0
			 * 0 0 [matrix                            ] 0 0
			 * 0 0 [boundary row for this element     ] 0 0
			 * 0 0 [interpolating row for this element] 0 0
			 * 0 0 interpolating row for below element  0 0
			 * 0 0 boundary row for below element       0 0
			 ************************************************************************/
			fourier_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
			fourier_solver (bases::master_solver <datatype> &i_solver, datatype& i_timestep, std::shared_ptr <bases::boundary <datatype>> i_boundary_0, std::shared_ptr <bases::boundary <datatype>> i_boundary_n);
			
			virtual ~fourier_solver () {}
			
			datatype *matrix_ptr () {
				return &horizontal_matrix [0];
			}
	
			void factorize ();
			void execute ();
		};
		
		template <class datatype>
		class laplace_solver : public bases::solver <datatype>
		{
		public:
			laplace_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
			laplace_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr);
			
			virtual ~laplace_solver () {}
			
			datatype *matrix_ptr () {
				return NULL;
			}
	
			void factorize ();
			void execute ();
		
		private:
			int n;
			int ldn;
			int m;
			datatype *data;
			datatype ex_pos_0, ex_pos_m;
			int flags;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;
			const datatype *pos_n, *pos_m;
			datatype *sub_ptr, *diag_ptr, *sup_ptr;
			int excess_0, excess_n, id, np;
			datatype *rhs_ptr;

			utils::messenger* messenger_ptr;
			
			std::vector <datatype> x;
			std::vector <datatype> sup, sub, diag, supsup; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <int> ipiv, xipiv;
		};
		
		template <class datatype>
		class incompressible_corrector : public bases::solver <datatype>
		{
		public:
			incompressible_corrector (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, datatype *i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags, int *i_component_x, int *i_component_z);
			incompressible_corrector (bases::master_solver <datatype> &i_solver, bases::master_solver <datatype> &i_solver_x, bases::master_solver <datatype> &i_solver_z, utils::messenger* i_messenger_ptr);
			
			virtual ~incompressible_corrector () {}
			
			datatype *matrix_ptr () {
				return NULL;
			}
	
			void factorize ();
			void execute ();
		
		private:
			int n;
			int ldn;
			int m;
			int ntop, nbot;
			datatype *data, *data_x, *data_z;
			datatype ex_pos_0, ex_pos_m, exx_pos_0, exx_pos_m, exxx_pos_0, exxx_pos_m;
			int flags;
			int *component_flags_x, *component_flags_z;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;
			const datatype *pos_n;
			datatype *pos_m;
			datatype *sub_ptr, *diag_ptr, *sup_ptr;
			int excess_0, excess_n, id, np;
			datatype *rhs_ptr;
			std::shared_ptr <bases::plan <datatype> > transform, transform_h, x_deriv, z_deriv;

			utils::messenger* messenger_ptr;
			
			std::vector <datatype> x, bufferl, bufferr;
			std::vector <datatype> data_temp, positions;
			std::vector <datatype> sup, sub, diag, supsup, matrix; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <int> ipiv, xipiv;
		};
	} /* fourier */
} /* two_d */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
