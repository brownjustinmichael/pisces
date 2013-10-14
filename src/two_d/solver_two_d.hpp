/*!**********************************************************************
 * \file solver_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-11.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_TWO_D_HPP_JW4BV4PS
#define SOLVER_TWO_D_HPP_JW4BV4PS

#include "../bases/solver.hpp"

namespace two_d
{
	enum solve_flags {
		x_solve = 0x20,
		y_solve = 0x40
	}
	
	namespace fourier
	{
		namespace chebyshev
		{
			class solver : public bases::solver
			{
			public:
				solver (bases::messenger* i_messenger_ptr, int i_n, int i_m, int i_excess_0, int i_excess_n, int i_n_iterations, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* positions, datatype *i_default_matrix, datatype *i_matrix, datatype* i_data_in, datatype* i_explicit_rhs, datatype* i_implicit_rhs, datatype* i_data_out = NULL, int i_flags = 0x00);
				
				virtual ~solver ();
			
			private:
				/* data */
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: SOLVER_TWO_D_HPP_JW4BV4PS */
