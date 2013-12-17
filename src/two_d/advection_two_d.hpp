/*!**********************************************************************
 * \file advection_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-16.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ADVECTION_TWO_D_HPP_GGR0NN1Q
#define ADVECTION_TWO_D_HPP_GGR0NN1Q

#include "plan_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace chebyshev
		{
			template 
			class advection : two_d::explicit_plan
			{
			public:
				advection (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype* i_vel_n, datatype *i_vel_m, datatype* i_data_in = NULL, datatype *i_data_in_2 = NULL) :
				explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, NULL),
				coeff (i_coeff) {}
				
				virtual ~advection () {}
				
				virtual void execute (int &element_flags) {
					for (int i = 0; i < m; ++i) {
						for (int j = 0; j < n; ++j) {
							
						}
					}
				}
			
			private:
				datatype coeff;
				using implicit_plan <datatype>::n;
				using implicit_plan <datatype>::m;
				using implicit_plan <datatype>::data_in;
				using implicit_plan <datatype>::data_out;
				using implicit_plan <datatype>::matrix_n;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ADVECTION_TWO_D_HPP_GGR0NN1Q */
