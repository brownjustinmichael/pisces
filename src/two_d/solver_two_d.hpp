/*!**********************************************************************
 * \file solver_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_TWO_D_HPP_JW4BV4PS
#define SOLVER_TWO_D_HPP_JW4BV4PS

#include "../bases/messenger.hpp"
#include "../bases/solver.hpp"
#include "element_two_d.hpp"
#include "plan_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		template <class datatype>
		class solver : public bases::solver <datatype>
		{
		public:
			solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, bases::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_aplha_n, datatype* i_data, int *i_element_flags, int *i_component_flags);
			
			virtual ~solver () {}
			
			virtual void reset () {
				utils::scale (ldn * m, 0.0, rhs_ptr (implicit_rhs));
				utils::scale (ldn * m, 0.0, rhs_ptr (explicit_rhs));
				utils::scale (ldn * m, 0.0, rhs_ptr (real_rhs));
			}
			
			datatype *matrix_ptr (int index = 0) {
				if (index == 0) {
					return &horizontal_matrix [0];
				} else {
					return &matrix [0];
				}
			}
	
			datatype *data_ptr () {
				return data;
			}
	
			datatype *rhs_ptr (int index = implicit_rhs) {
				if (index == implicit_rhs) {
					return &implicit_rhs_vec [0];
				} else if (index == explicit_rhs) {
					return &explicit_rhs_vec [0];
				} else if (index == real_rhs) {
					return &real_rhs_vec [0];
				} else {
					return NULL;
				}
			}
	
			bases::grid <datatype> *grid_ptr (int index = 0) {
				if (index == 0) {
					return &grid_n;
				} else {
					return &grid_m;
				}
			}
	
			using bases::solver <datatype>::element_flags;
			using bases::solver <datatype>::component_flags;
		
		private:
			void _factorize ();
			void _solve ();
			
			int n;
			int ldn;
			int m;
			datatype *data;
			int flags;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;

			bases::messenger* messenger_ptr;
	
			datatype& timestep; //!< A datatype reference to the current timestep
			datatype& alpha_0; //!< A datatype reference to the current edge_0 weight
			datatype& alpha_n; //!< A datatype reference to the current edge_n weight
			std::vector <datatype> values_0, values_n;

			datatype* positions;
			int temp_n;
			int excess_0; //!< The integer number of elements to recv from edge_0
			int excess_n; //!< The integer number of elements to recv from edge_n
			int ex_excess_0; //!< The integer number of elements to send to edge_0
			int ex_excess_n; //!< The integer number of elements to send to edge_n
			int ntop, nbot;

			datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component

			std::vector <datatype> explicit_rhs_vec;
			std::vector <datatype> implicit_rhs_vec;
			std::vector <datatype> real_rhs_vec;
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
			laplace_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, bases::messenger* i_messenger_ptr, datatype* i_data, int *i_element_flags, int *i_component_flags);
			
			virtual ~laplace_solver () {}
			
			virtual void reset () {
				utils::scale (ldn * m, 0.0, rhs_ptr (explicit_rhs));
				utils::scale (ldn * m, 0.0, rhs_ptr (real_rhs));
			}
			
			datatype *matrix_ptr (int index = 0) {
				return NULL;
			}
	
			datatype *data_ptr () {
				return data;
			}
	
			datatype *rhs_ptr (int index = explicit_rhs) {
				if (index == explicit_rhs) {
					return &explicit_rhs_vec [0];
				} else if (index == real_rhs) {
					return &real_rhs_vec [0];
				} else {
					return NULL;
				}
			}
	
			bases::grid <datatype> *grid_ptr (int index = 0) {
				if (index == 0) {
					return &grid_n;
				} else {
					return &grid_m;
				}
			}
	
		
		private:
			void _factorize ();
			void _solve ();
	
			int n;
			int ldn;
			int m;
			datatype *data;
			datatype ex_pos_0, ex_pos_m;
			int flags;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;
			datatype *pos_n, *pos_m, *sub_ptr, *diag_ptr, *sup_ptr;
			int excess_0, excess_n, id, np;

			bases::messenger* messenger_ptr;
			
			std::vector <datatype> x;
			std::vector <datatype> explicit_rhs_vec;
			std::vector <datatype> real_rhs_vec;
			std::vector <datatype> sup, sub, diag, supsup; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <int> ipiv, xipiv;
			
			std::shared_ptr <bases::plan <datatype> > transform;
			
		};
		
		template <class datatype>
		class divergence_solver : public bases::solver <datatype>
		{
		public:
			divergence_solver (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_x, datatype *i_data_z, int *i_element_flags, int *i_component_flags);
			
			virtual ~divergence_solver () {}
			
			virtual void reset () {}
			
			datatype *matrix_ptr (int index = 0) {
				return NULL;
			}
	
			datatype *data_ptr () {
				return data_x;
			}
	
			datatype *rhs_ptr (int index = explicit_rhs) {
				return NULL;
			}
	
			bases::grid <datatype> *grid_ptr (int index = 0) {
				if (index == 0) {
					return &grid_n;
				} else {
					return &grid_m;
				}
			}
	
		
		private:
			void _factorize ();
			void _solve ();
	
			int n;
			int ldn;
			int m;
			datatype *data_x, *data_z, *pos_m, scalar;
			int flags;
			bases::grid <datatype> &grid_n;
			bases::grid <datatype> &grid_m;
		};
	} /* fourier */
} /* two_d */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
