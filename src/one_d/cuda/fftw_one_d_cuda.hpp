/*!***********************************************************************
 * \file fftw_one_d_cuda.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFTW_ONE_D_CUDA_HPP_G5118SR0
#define FFTW_ONE_D_CUDA_HPP_G5118SR0

#include "../../config.hpp"
#include "../../bases/plan.hpp"

namespace one_d
{
	namespace cuda
	{
		/*!*******************************************************************
		 * An implementation of the transform class using FFTW3.
		 * 
		 * \brief \copybrief bases::explicit_plan <double>
		 *********************************************************************/
		class fftw_cosine : public bases::plan <double>
		{
		public:
			/*!*******************************************************************
			 * \copydoc bases::explicit_plan <double>::explicit_plan ()
			 *********************************************************************/
			fftw_cosine (bases::element <double>* i_element_ptr, int i_n, double* i_data_dev);
		
			virtual ~fftw_cosine ();
		
			/*!*******************************************************************
			 * \copydoc bases::explicit_plan <double>::execute ()
			 *********************************************************************/
			void execute ();
	
		private:
			int n;
			double scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
			void* data_real;
			void* data_complex;
			void* cu_plan;
			// fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		};
		
		class transfer : public bases::plan <double>
		{
		public:
			transfer (bases::element <double>* i_element_ptr, int i_n, double* i_data_dev, double* i_data);
			
			virtual ~transfer () {}
			
			virtual void execute ();
		
		private:
			int n;
			double* data_dev;
			double* data;
		};
	} /* cuda */
} /* one_d */

#endif /* end of include guard: FFTW_ONE_D_CUDA_HPP_G5118SR0 */
