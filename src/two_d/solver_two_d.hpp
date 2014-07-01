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
#include "plan_two_d.hpp"

namespace two_d
{
	template <class datatype>
	class master_solver : public bases::master_solver <datatype>
	{
	private:
		int n;
		int ldn;
		int m;
		bases::grid <datatype> &grid_n;
		bases::grid <datatype> &grid_m;
		
		std::shared_ptr <bases::solver <datatype> > x_solver;
		std::shared_ptr <bases::solver <datatype> > z_solver;
		
		std::vector <datatype> explicit_rhs_vec;
		std::vector <datatype> implicit_rhs_vec;
		std::vector <datatype> real_rhs_vec;
		
		datatype *explicit_rhs_ptr = NULL;
		datatype *implicit_rhs_ptr = NULL;
		datatype *real_rhs_ptr = NULL;
		
		using bases::master_solver <datatype>::data;
		using bases::master_solver <datatype>::element_flags;
		using bases::master_solver <datatype>::component_flags;
		
	public:
		master_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_data, int *i_element_flags, int *i_component_flags) : bases::master_solver <datatype> (i_data, i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), grid_n (i_grid_n), grid_m (i_grid_m) {
			// *component_flags |= z_solve;
		}
		
		virtual ~master_solver () {}
		
		bases::grid <datatype> *grid_ptr (int index = 0) {
			if (index == 0) {
				return &grid_n;
			} else {
				return &grid_m;
			}
		}
		
		datatype *rhs_ptr (int index = implicit_rhs) {
			if (index == implicit_rhs) {
				if (!implicit_rhs_ptr) {
					implicit_rhs_vec.resize (ldn * m);
					implicit_rhs_ptr = &implicit_rhs_vec [0];
				}
				return implicit_rhs_ptr;
			} else if (index == explicit_rhs) {
				if (!explicit_rhs_ptr) {
					explicit_rhs_vec.resize (ldn * m);
					explicit_rhs_ptr = &explicit_rhs_vec [0];
				}
				return explicit_rhs_ptr;
			} else if (index == real_rhs) {
				if (!real_rhs_ptr) {
					real_rhs_vec.resize (ldn * m);
					real_rhs_ptr = &real_rhs_vec [0];
				}
				return real_rhs_ptr;
			} else {
				return NULL;
			}
		}
		
		datatype *matrix_ptr (int index = 0) {
			if (index == 0) {
				DEBUG ("X POINTER IS " << x_solver->matrix_ptr (index));
				return x_solver->matrix_ptr (index);
			} else {
				DEBUG ("Z POINTER IS " << z_solver->matrix_ptr (index));
				return z_solver->matrix_ptr (index);
			}
		}
		
		virtual void reset () {
			if (implicit_rhs_ptr) {
				utils::scale (ldn * m, 0.0, implicit_rhs_ptr);
			}
			if (explicit_rhs_ptr) {
				utils::scale (ldn * m, 0.0, explicit_rhs_ptr);
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
			if (!(flags & not_x_solver)) {
				x_solver = i_solver;
			}
			if (!(flags & not_z_solver)) {
				z_solver = i_solver;
			}
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
	
	template <class datatype>
	class boundary
	{
	public:
		boundary () {}
		
		virtual ~boundary () {}
		
		virtual void calculate_rhs () = 0;
	};
	
	template <class datatype>
	class communicating_boundary : public boundary <datatype>
	{
	private:
		datatype alpha;
		utils::messenger *messenger_ptr;
		// std::vector <datatype> boundary_positions;
		int n_boundary_in;
		int n_boundary_out;
		int out_id;
		
	public:
		communicating_boundary (utils::messenger *i_messenger_ptr, int n, int &i_n_boundary_in, int &i_n_boundary_out, const datatype *i_positions, std::vector <datatype> &boundary_positions, bool top, int &bound) : messenger_ptr (i_messenger_ptr) {
			n_boundary_in = i_n_boundary_in;
			out_id = messenger_ptr->get_id () + (top ? 1 : -1);
			messenger_ptr->template send <int> (1, &i_n_boundary_in, out_id, 0);
			messenger_ptr->template recv <int> (1, &i_n_boundary_out, out_id, 0);
			n_boundary_out = i_n_boundary_out;
			boundary_positions.resize (i_n_boundary_out);
			messenger_ptr->template send <datatype> (i_n_boundary_in, i_positions, out_id, 0);
			messenger_ptr->template recv <datatype> (i_n_boundary_out, &boundary_positions [0], out_id, 0);
			bound = 1;
		}
		
		virtual ~communicating_boundary () {}
		
		virtual void calculate_rhs () {
			
		}
	};
	
	namespace fourier
	{
		template <class datatype>
		class collocation_solver : public bases::solver <datatype>
		{
		public:
			collocation_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype *i_implicit_rhs, datatype *i_explicit_rhs, datatype *i_real_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
			collocation_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_aplha_n);
			
			virtual ~collocation_solver () {}
			
			datatype *matrix_ptr (int index = 0) {
				if (index == 0) {
					return &matrix [0];
				} else {
					return &horizontal_matrix [0];
				}
			}
	
			void factorize ();
			void execute ();
			
			using bases::solver <datatype>::element_flags;
			using bases::solver <datatype>::component_flags;
		
		private:
			int n;
			int ldn;
			int m;
			datatype *data;
			int flags;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;

			utils::messenger* messenger_ptr;
	
			datatype& timestep; //!< A datatype reference to the current timestep
			datatype& alpha_0; //!< A datatype reference to the current edge_0 weight
			datatype& alpha_n; //!< A datatype reference to the current edge_n weight
			std::vector <datatype> values_0, values_n;
			datatype *implicit_rhs_vec, *explicit_rhs_vec, *real_rhs_vec;

			const datatype* positions;
			int temp_n;
			int excess_0; //!< The integer number of elements to recv from edge_0
			int excess_n; //!< The integer number of elements to recv from edge_n
			int ex_excess_0; //!< The integer number of elements to send to edge_0
			int ex_excess_n; //!< The integer number of elements to send to edge_n
			int ntop, nbot;

			datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component

			std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <datatype> positions_0; //!< A datatype vector of excess positions from edge_0
			std::vector <datatype> positions_n; //!< A datatype vector of excess positions from edge_n
			std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
			std::vector <datatype> boundary_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
			std::vector <datatype> previous_rhs;
			std::vector <int> ns;
			std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
			std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
			std::vector <datatype> matrix;
			std::vector <datatype> horizontal_matrix;
			std::vector <datatype> factorized_horizontal_matrix;
			
			std::shared_ptr <bases::plan <datatype> > transform;
		};
		
		template <class datatype>
		class laplace_solver : public bases::solver <datatype>
		{
		public:
			laplace_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, utils::messenger* i_messenger_ptr, datatype* i_implicit_rhs, datatype *i_explicit_rhs, datatype *i_real_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
			laplace_solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr);
			
			virtual ~laplace_solver () {}
			
			datatype *matrix_ptr (int index = 0) {
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
			datatype *implicit_rhs_vec, *explicit_rhs_vec, *real_rhs_vec;

			utils::messenger* messenger_ptr;
			
			std::vector <datatype> x;
			std::vector <datatype> sup, sub, diag, supsup; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <int> ipiv, xipiv;
			
			std::shared_ptr <bases::plan <datatype> > transform;
			
		};
		
		template <class datatype>
		class divergence_solver : public bases::solver <datatype>
		{
		public:
			divergence_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags);
			divergence_solver (bases::master_solver <datatype> &i_solver, datatype *i_data_z);
			
			virtual ~divergence_solver () {}
			
			datatype *matrix_ptr (int index = 0) {
				return NULL;
			}
			
			void factorize ();
			void execute ();
			
		private:
			int n;
			int ldn;
			int m;
			const datatype *pos_m;
			datatype *data_x, *data_z, scalar;
			int flags;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;
		};
	} /* fourier */
} /* two_d */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
