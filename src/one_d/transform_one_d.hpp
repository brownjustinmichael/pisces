/*!***********************************************************************
 * \file transform_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFTW_HPP_P3TP70YE
#define FFTW_HPP_P3TP70YE

#ifdef _CUDA
#include "cuda/transform_one_d_cuda.hpp"
#endif /* _CUDA */
#include <fftw3.h>
#include <memory>
#include "logger/logger.hpp"
#include "../bases/transform.hpp"
#include "plan_one_d.hpp"
#include "linalg/utils.hpp"

namespace one_d
{
	/*!*******************************************************************
	 * An implementation of the transform class using FFTW3.
	 * 
	 * \brief \copybrief bases::transform
	 *********************************************************************/
	template <class datatype>
	class transform : public bases::plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::transform::transform ()
		 *********************************************************************/
		transform (bases::grid <datatype> &i_grid, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags);
	
		virtual ~transform () {
			// printf ("Destroying one_d transform\n");
		}
	
		/*!*******************************************************************
		 * \copydoc bases::transform::execute ()
		 *********************************************************************/
		void execute ();
	
	protected:
		int n;
		datatype *data_in;
		datatype *data_out;
	
		datatype scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
		fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		fftwf_plan fourier_plan_float; //!< The fftw_plan object to be executed
	};
	
	template <class datatype>
	class master_transform : public bases::master_transform <datatype>
	{
	public:
		master_transform (bases::grid <datatype> &i_grid, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags):
		bases::master_transform <datatype> (element_flags, component_flags),
		ldn (i_grid.get_ld ()),
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in) {
			data.resize (ldn);
			if (i_flags & forward_vertical) {
				forward_transform = std::shared_ptr <bases::plan <datatype>> (new one_d::transform <datatype> (i_grid, &data [0], &data [0], 0x00, element_flags, component_flags));
			}
			if (i_flags & inverse_vertical) {
				if (forward_transform) {
					inverse_transform = forward_transform;
				} else {
					inverse_transform = std::shared_ptr <bases::plan <datatype>> (new one_d::transform <datatype> (i_grid, &data [0], &data [0], inverse, element_flags, component_flags));
				}
			}
		}
		
		void _transform (int flags) {
			if (flags & forward_vertical) {
				if (!(*component_flags & transformed_vertical) && forward_transform) {
					forward_transform->execute ();
				} else if (*component_flags & transformed_vertical) {
					WARN ("Forward transform attempted on spectral data: ignoring.")
				}
			}
			if (flags & inverse_vertical) {
				if ((*component_flags & transformed_vertical) && inverse_transform) {
					inverse_transform->execute ();
				} else if (!(*component_flags & transformed_vertical)) {
					WARN ("Inverse transform attempted on spectral data: ignoring.")
				}
			}
		}
		
		void write () {
			utils::copy (ldn, data_in, &data [0]);
		}

		void read () {
			utils::copy (ldn, &data [0], data_out);
		}
	private:
		int ldn;
		datatype *data_in, *data_out;
		std::vector <datatype> data;
		
		std::shared_ptr <bases::plan <datatype>> forward_transform;
		std::shared_ptr <bases::plan <datatype>> inverse_transform;
		
		using bases::master_transform <datatype>::component_flags;
	};
} /* one_d */

#endif /* end of include guard: FFTW_HPP_P3TP70YE */
