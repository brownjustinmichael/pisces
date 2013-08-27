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

namespace cuda
{
	namespace one_d
	{
		/*!*******************************************************************
		 * An implementation of the transform class using FFTW3.
		 * 
		 * \brief \copybrief bases::explicit_plan <double>
		 *********************************************************************/
		template <class datatype>
		class fftw_cosine : public bases::plan <datatype>
		{
		public:
			/*!*******************************************************************
			 * \copydoc bases::explicit_plan <double>::explicit_plan ()
			 *********************************************************************/
			fftw_cosine (bases::element <datatype>* i_element_ptr, int i_n, datatype* i_data_dev);
		
			virtual ~fftw_cosine ();
		
			/*!*******************************************************************
			 * \copydoc bases::explicit_plan <double>::execute ()
			 *********************************************************************/
			void execute ();
	
		private:
			using bases::plan <datatype>::flags_ptr;

			int n;
			datatype scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
			void* data_real;
			void* data_complex;
			void* cu_plan;
			// fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		};
		
		template <class datatype>
		class transfer : public bases::plan <datatype>
		{
		public:
			transfer (bases::element <datatype>* i_element_ptr, int i_n, datatype* i_data_dev, datatype* i_data);
			
			virtual ~transfer () {}
			
			virtual void execute ();
		
		private:
			using bases::plan <datatype>::element_ptr;
			using bases::plan <datatype>::flags_ptr;
			
			int n;
			datatype* data_dev;
			datatype* data;
		};
	} /* one_d */
} /* cuda */

#endif /* end of include guard: FFTW_ONE_D_CUDA_HPP_G5118SR0 */
