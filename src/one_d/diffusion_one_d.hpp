/*!***********************************************************************
 * \file diffusion_one_d.hpp
 * Spectral Element
 * 
 * This file provides several implementations of implicit and explicit 
 * methods for solving diffusion.
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_ONE_D_H_1AIIK1RA
#define DIFFUSION_ONE_D_H_1AIIK1RA

#include <memory>
#include <vector>
#include "plan_one_d.hpp"
#include "../bases/grid.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	template <class datatype>
	class diffusion : public implicit_plan <datatype>
	{
	public:
		diffusion (bases::grid <datatype> &i_grid, datatype i_coeff, datatype i_alpha, datatype *i_matrix, datatype* i_data_in, datatype* i_data_out = NULL);
		
		virtual ~diffusion () {
			// printf ("Destroying one_d diffusion\n");
		}
		
		void execute (bases::flags &element_flags);
	
	private:
		datatype coeff;
		datatype alpha;
		using implicit_plan <datatype>::n;
		using implicit_plan <datatype>::data_in;
		using implicit_plan <datatype>::data_out;
		using implicit_plan <datatype>::matrix;
		using implicit_plan <datatype>::grid;
	};
	
	template <class datatype>
	class nonlinear_diffusion : public explicit_plan <datatype>
	{
	public:
		nonlinear_diffusion (bases::grid <datatype> &i_grid, datatype i_coeff, datatype* i_data_in, datatype* i_data_out = NULL);
		
		virtual ~nonlinear_diffusion () {
			// printf ("Destroying one_d nonlinear_diffusion\n");
		}
		
		void execute (bases::flags &element_flags);
	
	private:
		datatype coeff;
		datatype alpha;
		using explicit_plan <datatype>::n;
		using explicit_plan <datatype>::data_in;
		using explicit_plan <datatype>::data_out;
		using explicit_plan <datatype>::grid;
	};
} /* oned */

#endif /* end of include guard: DIFFUSION_ONE_D_H_1AIIK1RA */
