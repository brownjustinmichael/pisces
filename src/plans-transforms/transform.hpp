/*!**********************************************************************
 * \file transform_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_TWO_D_HPP_RAXYBFTC
#define TRANSFORM_TWO_D_HPP_RAXYBFTC

#include "logger/logger.hpp"
#include "plans/plan.hpp"
#include "transform.hpp"
#include <fftw3.h>
#include "linalg/utils.hpp"

namespace plans
{
	// class fftw_configure
	// {
	// public:
	// 	fftw_configure ();
	// 	
	// 	virtual ~fftw_configure ();
	// 	
	// private:
	// 	bool threads;
	// };
	
	template <class datatype>
	class horizontal_transform : public plans::plan <datatype>
	{
	public:
		/*!**********************************************************************
		 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
		 ************************************************************************/
		horizontal_transform (int n, int m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);

		horizontal_transform (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);
		
		virtual ~horizontal_transform () {}
		
		virtual void execute ();
	
	protected:
		int n;
		int m;
		datatype *data_in;
		datatype *data_out;
		
		int flags;
		int threads;
		datatype scalar;
		std::vector <fftw_plan> plans;
		std::vector <fftwf_plan> plans_float;
		fftw_iodim major_iodim;
		fftw_iodim iodim;
	};
	
	template <class datatype>
	class vertical_transform : public plans::plan <datatype>
	{
	public:
		/*!**********************************************************************
		 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
		 ************************************************************************/
		vertical_transform (int n, int m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);

		vertical_transform (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);
		
		virtual ~vertical_transform () {}
		
		virtual void execute ();
	
	protected:
		int n;
		int m;
		datatype *data_in;
		datatype *data_out;
		
		int flags;
		int threads;
		datatype scalar;
		std::vector <fftw_plan> plans;
		std::vector <fftwf_plan> plans_float;
	};
} /* plans */

#endif /* end of include guard: TRANSFORM_TWO_D_HPP_RAXYBFTC */
