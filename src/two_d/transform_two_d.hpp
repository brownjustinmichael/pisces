/*!**********************************************************************
 * \file transform_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_TWO_D_HPP_RAXYBFTC
#define TRANSFORM_TWO_D_HPP_RAXYBFTC

#include "../config.hpp"
#include "plan_two_d.hpp"
#include <fftw3.h>

namespace two_d
{
	enum transform_flags {
		ignore_m = 0x01,
		inverse = 0x02
	};
	
	
	/*
		TODO 
	*/
	// class fftw_configure
	// {
	// public:
	// 	fftw_configure (int n_threads);
	// 	
	// 	virtual ~fftw_configure ();
	// 
	// private:
	// 	int n_threads;
	// };
	
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class horizontal_transform : public bases::plan <datatype>
			{
			public:
				/*!**********************************************************************
				 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
				 ************************************************************************/
				horizontal_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out = NULL, int i_flags = 0x00);
				
				virtual ~horizontal_transform () {}
				
				virtual void execute (int &element_flags, int &component_flags);
			
			protected:
				int n;
				int m;
				datatype *data_in;
				datatype *data_out;
				
				int flags;
				datatype scalar;
				fftw_plan x_plan;
				fftwf_plan x_plan_float;
				fftw_iodim major_iodim;
				fftw_iodim iodim;
			};
			
			template <class datatype>
			class vertical_transform : public bases::plan <datatype>
			{
			public:
				/*!**********************************************************************
				 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
				 ************************************************************************/
				vertical_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out = NULL, int i_flags = 0x00);
				
				virtual ~vertical_transform () {}
				
				virtual void execute (int &element_flags, int &component_flags);
			
			protected:
				int n;
				int m;
				datatype *data_in;
				datatype *data_out;
				
				int flags;
				datatype scalar;
				fftw_plan z_plan;
				fftwf_plan z_plan_float;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: TRANSFORM_TWO_D_HPP_RAXYBFTC */
