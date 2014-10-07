/*!**********************************************************************
 * \file fourier_solver.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FOURIER_SOLVER_HPP_6DE77A57
#define FOURIER_SOLVER_HPP_6DE77A57

#include "../solver.hpp"
#include "../boundary.hpp"

namespace plans
{
	template <class datatype>
	class fourier_solver : public plans::solver <datatype>
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
		
		std::shared_ptr <plans::boundary <datatype>> boundary_0, boundary_n;
		
		int inner_m;
		int ex_overlap_0;
		int overlap_0;
		int ex_overlap_n;
		int overlap_n;
		int lda;
		const datatype *pos_m;
		
		using plans::solver <datatype>::element_flags;
		using plans::solver <datatype>::component_flags;
		
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
		fourier_solver (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype& i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
				
		virtual ~fourier_solver () {}
		
		datatype *matrix_ptr () {
			return &horizontal_matrix [0];
		}

		void factorize ();
		void execute ();
		
		class factory : public plans::solver <datatype>::factory
		{
		private:
			datatype &timestep;
			std::shared_ptr <plans::boundary <datatype>> boundary_0, boundary_n;

		public:
			factory (datatype &i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n) : timestep (i_timestep), boundary_0 (i_boundary_0), boundary_n (i_boundary_n) {}
			virtual ~factory () {}
			
			virtual std::shared_ptr <plans::solver <datatype>> instance (plans::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				return std::shared_ptr <plans::solver <datatype>> (new fourier_solver (*grids [0], *grids [1], timestep, boundary_0, boundary_n, i_rhs, i_data, i_element_flags, i_component_flags));
			}
		};
	};
} /* plans */

#endif /* end of include guard: FOURIER_SOLVER_HPP_6DE77A57 */