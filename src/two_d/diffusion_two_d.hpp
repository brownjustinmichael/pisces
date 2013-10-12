/*!**********************************************************************
 * \file diffusion_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_TWO_D_HPP_YHECX9VS
#define DIFFUSION_TWO_D_HPP_YHECX9VS

#include "../bases/plan.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class diffusion : public implicit_plan
			{
			public:
				diffusion (int i_n, int i_m, datatype i_coeff, datatype i_alpha, bases::grid <datatype>* i_grid_n, bases::grid <datatype>* i_grid_m, datatype* i_data_in, datatype* i_matrix, datatype* i_data_out = NULL);
				
				virtual ~diffusion ();
			
			private:
				/* data */
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
