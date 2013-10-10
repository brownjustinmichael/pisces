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
#include "../bases/plan.hpp"
#include "../bases/collocation.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	template <class datatype>
	class diffusion : public bases::implicit_plan <datatype>
	{
	public:
		diffusion (int i_n, datatype i_coeff, datatype i_alpha, bases::collocation_grid <datatype>* i_grid, datatype* i_data_in, datatype* i_matrix, datatype* i_data_out = NULL, int i_flags = 0x0);
		
		virtual ~diffusion () {}
		
		void execute ();
	
	private:
		datatype coeff;
		datatype alpha;
		int flags;
		using bases::implicit_plan <datatype>::n;
		using bases::implicit_plan <datatype>::data_in;
		using bases::implicit_plan <datatype>::data_out;
		using bases::implicit_plan <datatype>::matrix;
		using bases::implicit_plan <datatype>::grid;
	};
} /* oned */

#endif /* end of include guard: DIFFUSION_ONE_D_H_1AIIK1RA */
