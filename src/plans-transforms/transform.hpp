/*!**********************************************************************
 * \file transform.hpp
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

/*!**********************************************************************
 * \namespace plans::transforms
 * 
 * \brief A namespace that contains the transforms and transform controllers
 ************************************************************************/
namespace plans
{
	namespace transforms
	{
		/*!**********************************************************************
		 * \brief A plan that transforms the data horizontally via Fourier transform
		 ************************************************************************/
		template <class datatype>
		class horizontal : public plans::plan <datatype>
		{
		protected:
			using plans::plan <datatype>::data_in;
			using plans::plan <datatype>::data_out;

			int n; //!< The number of horizontal data points
			int m; //!< The number of vertical data points
			int flags; //!< A set of flags for setting up the transform
			int threads; //!< The integer number of threads to use
			
			datatype scalar; //!< The scalar by which to scale the data after transform
			std::vector <fftw_plan> plans; //!< A vector of fftw_plan objects
			std::vector <fftwf_plan> plans_float; //!< A vector of fftwf_plan objects (for single precision)
			fftw_iodim major_iodim; //!< A dimensional object needed by fftw
			fftw_iodim iodim; //!< A dimensional object needed by fftw
			
		private:
			/**
			 * @brief Initialize the transform plan
			 * @details This member initializes the transform plan, including the dimensionalities and the FFTW plans that work beneath the hood.
			 * 
			 * @param i_n The number of horizontal data points
			 * @param i_m The number of vertical data points
			 * @param i_data_in The pointer to the input data
			 * @param i_data_out The pointer to the output data
			 * @param i_flags A set of flags for setting up the transform (e.g. inverse)
			 * @param i_threads The number of threads over which to split the transform
			 */
			void init (int i_n, int i_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int i_threads = 0);
			
		public:
			/*!**********************************************************************
			 * \param i_data_in A reference to the variable object of the input data
			 * \param i_data_out A reference to the variable object of the output data
			 * \param i_flags A set of flags for setting up the transform (e.g. inverse)
			 * \param i_threads The number of threads over which to split the transform
			 * @param state_in The integer input state of the data (e.g. real_real, real_spectral)
			 * @param state_out The integer output state of the data (e.g. real_real, real_spectral)
			 * 
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			horizontal (grids::variable &i_data_in, grids::variable &i_data_out, int state_in = 0, int state_out = 0, int i_flags = 0x00, int i_threads = 0) : 
			plans::plan <datatype> (i_data_in, i_data_out, state_in, state_out) {
				init (i_data_in.get_grid (0).get_n (), i_data_in.get_grid (1).get_n (), data_in, data_out, i_flags, i_threads);
			}

			/*!**********************************************************************
			 * @brief A shorthand constructor for an in-place transform
			 * 
			 * \param i_data_in A reference to the variable object of the input data
			 * \param i_flags A set of flags for setting up the transform (e.g. inverse)
			 * \param i_threads The number of threads over which to split the transform
			 * @param state_in The integer input state of the data (e.g. real_real, real_spectral)
			 * @param state_out The integer output state of the data (e.g. real_real, real_spectral)
			 * 
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			horizontal (grids::variable &i_data_in, int state_in = 0, int state_out = -1, int i_flags = 0x00, int i_threads = 0) : plans::plan <datatype> (i_data_in, i_data_in, state_in, state_out >= 0 ? state_out : state_in) {
				init (i_data_in.get_grid (0).get_n (), i_data_in.get_grid (1).get_n (), data_in, data_out, i_flags, i_threads);
			}
			
			virtual ~horizontal () {}
			
			/*!**********************************************************************
			 * \copydoc plan::execute
			 ************************************************************************/
			virtual void execute ();
		};
		
		/*!**********************************************************************
		 * \brief A plan that transforms the data vertically via Chebyshev/Cosine transform
		 ************************************************************************/
		template <class datatype>
		class vertical : public plans::plan <datatype>
		{
		protected:
			using plans::plan <datatype>::data_in;
			using plans::plan <datatype>::data_out;

			int n; //!< The number of data points in the horizontal direction
			int m; //!< The number of data points in the vertical direction
			int flags; //!< A set of flags for setting up the transform
			int threads; //!< The integer number of threads to use
			
			datatype scalar; //!< The scalar by which to scale the data after transform
			std::vector <fftw_plan> plans; //!< A vector of fftw_plan objects
			std::vector <fftwf_plan> plans_float; //!< A vector of fftwf_plan objects (for single precision)
			
		private:
			void init (int i_n, int i_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int i_threads = 0);
			
		public:
			/*!**********************************************************************
			 * \param i_data_in A reference to the variable object of the input data
			 * \param i_data_out A reference to the variable object of the output data
			 * \param i_flags A set of flags for setting up the transform (e.g. inverse)
			 * \param i_threads The number of threads over which to split the transform
			 * @param state_in The integer input state of the data (e.g. real_real, real_spectral)
			 * @param state_out The integer output state of the data (e.g. real_real, real_spectral)
			 * 
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			vertical (grids::variable &i_data_in, grids::variable &i_data_out, int state_in = 0, int state_out = 0, int i_flags = 0x00, int i_threads = 0) : 
			plans::plan <datatype> (i_data_in, i_data_out, state_in, state_out) {
				init (i_data_in.get_grid (0).get_n (), i_data_in.get_grid (1).get_n (), data_in, data_out, i_flags, i_threads);
			}

			/*!**********************************************************************
			 * @brief A shorthand constructor to use for an in-place transform
			 * 
			 * \param i_data_in A reference to the variable object of the input data
			 * \param i_flags A set of flags for setting up the transform (e.g. inverse)
			 * \param i_threads The number of threads over which to split the transform
			 * @param state_in The integer input state of the data (e.g. real_real, real_spectral)
			 * @param state_out The integer output state of the data (e.g. real_real, real_spectral)
			 * 
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)			 ************************************************************************/
			vertical (grids::variable &i_data_in, int state_in = 0, int state_out = -1, int i_flags = 0x00, int i_threads = 0) : plans::plan <datatype> (i_data_in, i_data_in, state_in, state_out >= 0 ? state_out : state_in) {
				init (i_data_in.get_grid (0).get_n (), i_data_in.get_grid (1).get_n (), data_in, data_out, i_flags, i_threads);
			}
			
			virtual ~vertical () {}
			
			/*!**********************************************************************
			 * \copydoc plan::execute
			 ************************************************************************/
			virtual void execute ();
		};
	} /* transforms */
} /* plans */

#endif /* end of include guard: TRANSFORM_TWO_D_HPP_RAXYBFTC */
