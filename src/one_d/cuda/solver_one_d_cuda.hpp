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

namespace cuda
{
	namespace one_d
	{
		template <class datatype>
		class solver : public ::one_d::solver <datatype>
		{
		public:
			solver (bases::messenger* i_messenger_ptr, int i_n, int i_excess_0, int i_excess_n, datatype& i_timestep, datatype& i_alpha_0, datatype& i_alpha_n, datatype* i_positions, datatype *i_default_matrix, datatype *i_matrix, datatype* i_data_in, datatype* i_rhs, datatype* i_data_out = NULL, int i_flags = 0x00);
			
			virtual ~solver ();
			
			virtual void execute ();
			
		protected:
			virtual void _factorize ();
			
			using ::one_d::solver <datatype>::n;
			using ::one_d::solver <datatype>::ipiv;
			using ::one_d::solver <datatype>::data_in;
			using ::one_d::solver <datatype>::data_out;
			using ::one_d::solver <datatype>::rhs;
			using ::one_d::solver <datatype>::alpha_0;
			using ::one_d::solver <datatype>::alpha_n;
			using ::one_d::solver <datatype>::timestep;
			using ::one_d::solver <datatype>::factorized_matrix;
			
			int* ipiv_dev;
			datatype* factorized_matrix_dev;
		};
	} /* one_d */
} /* cuda */

#endif /* end of include guard: SOLVER_ONE_D_HPP_YNU1VLWZ */
