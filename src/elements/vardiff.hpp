/*!**********************************************************************
 * \file vardiff.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VARDIFF_TWO_D_HPP_OQ800X4X
#define VARDIFF_TWO_D_HPP_OQ800X4X

#include "mpi/messenger.hpp"

#include "boussinesq.hpp"

namespace pisces
{
	/*!**********************************************************************
	 * \brief An element built to use the Boussinesq approximation
	 ************************************************************************/
	class vardiff_element : public boussinesq_element
	{
	private:
		std::map <std::string, std::vector <double>> diffusion; //!< A vector of diffusion data, for background diffusion
		using boussinesq_element::equations;
		using boussinesq_element::ptr;
		using boussinesq_element::m;
		using boussinesq_element::messenger_ptr;
		
	public:
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		vardiff_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~vardiff_element () {}

		static std::shared_ptr <element> instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) {
			return std::shared_ptr <element> (new vardiff_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags));
		}
		
		/**
		 * @brief An alternate rezone merit function that seeks the regions with strongest derivatives
		 * 
		 * @param element_ptr A pointer to the element in question
		 * @param virt A virtual file on which to operate the merit function
		 * 
		 * @return A value to be minimized in order to get the boundaries to occur at strongest derivatives
		 */
		static double rezone_merit (element *element_ptr, formats::virtual_file *virt) {
			double *temp = &(virt->index <double> ("temperature"));
			double *vel = &(virt->index <double> ("z_velocity"));
			double *pos = &(virt->index <double> ("z"));
			double value = 0.0;
			for (int j = 0; j < virt->dims ["z"] [1] - 1; ++j) {
				value += -(vel [j] - vel [j + 1]) * (vel [j] - vel [j + 1]) - 1.0e3 * (temp [j] - temp [j + 1]) * (temp [j] - temp [j + 1]) + 1.0e-5 / ((pos [j] - pos [j + 1]) * (pos [j] - pos [j + 1]));
			}
			return value;
		}
	};
} /* pisces */

#endif /* end of include guard: VARDIFF_TWO_D_HPP_OQ800X4X */
