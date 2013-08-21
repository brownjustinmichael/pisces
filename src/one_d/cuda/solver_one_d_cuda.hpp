/*!**********************************************************************
 * \file solver_one_d.hpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_ONE_D_HPP_YNU1VLWZ
#define SOLVER_ONE_D_HPP_YNU1VLWZ

#include "../solver_one_d.hpp"
#include "element_one_d_cuda.hpp"

namespace one_d
{
	namespace cuda
	{
		template <class datatype>
		class solver : one_d::solver <datatype>
		{
		public:
			solver (bases::element <datatype>* i_element_ptr, int i_n, int i_excess_0, int i_excess_n, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype *i_default_matrix, datatype *i_matrix, int i_name_in, int i_name_rhs, int i_name_out = null, int i_flags = 0x00);
			
			virtual ~solver () {}
			
			virtual void execute ();
			
		protected:
			virtual void _factorize ();
			
			datatype* factorized_matrix;
			datatype* factorized_matrix_dev;
		};
	} /* cuda */
} /* one_d */

#endif /* end of include guard: SOLVER_ONE_D_HPP_YNU1VLWZ */
