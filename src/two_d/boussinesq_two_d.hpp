/*!**********************************************************************
 * \file boussinesq_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUSSINESQ_TWO_D_HPP_OQ800X4X
#define BOUSSINESQ_TWO_D_HPP_OQ800X4X

#include "element_two_d.hpp"

namespace two_d
{
	namespace fourier
	{
		namespace cosine
		{	
			template <class datatype>
			class boussinesq_element : public element <datatype>
			{
			public:
				boussinesq_element (bases::axis i_axis_n, bases::axis i_axis_m, int i_name, io::parameters& i_params, bases::messenger* i_messenger_ptr, int i_element_flags);
				
				virtual ~boussinesq_element () {}
			
				datatype calculate_timestep (int i, int j, io::virtual_dump *dump = NULL);
				
				using element <datatype>::ptr;
			private:
				using element <datatype>::element_flags;
				using element <datatype>::params;
				using element <datatype>::initialize;
				using element <datatype>::n;
				using element <datatype>::m;
				using element <datatype>::name;
				using element <datatype>::normal_stream;
				using element <datatype>::transform_stream;
				using element <datatype>::cell_n;
				using element <datatype>::cell_m;
				using element <datatype>::grids;
				using element <datatype>::matrix_ptr;
				using element <datatype>::messenger_ptr;
				using element <datatype>::timestep;
				using element <datatype>::duration;
				using element <datatype>::alpha_0;
				using element <datatype>::alpha_n;
				using element <datatype>::solvers;
				using element <datatype>::scalars;
				using element <datatype>::scalar_names;
				
				datatype advection_coeff, cfl, *x_ptr, *z_ptr, *x_vel_ptr, *z_vel_ptr;
			};
		} /* cosine */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: BOUSSINESQ_TWO_D_HPP_OQ800X4X */
