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
	template <class datatype>
	class vardiff_element : public boussinesq_element <datatype>
	{
	private:
		std::map <std::string, std::vector <datatype>> diffusion; //!< A vector of diffusion data, for background diffusion
		using boussinesq_element <datatype>::equations;
		using boussinesq_element <datatype>::ptr;
		using boussinesq_element <datatype>::m;
		using boussinesq_element <datatype>::messenger_ptr;
		
	public:
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		vardiff_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~vardiff_element () {}
		
		static double rezone_merit (element <datatype> *element_ptr, formats::virtual_file *virt) {
			datatype *vel = &(virt->index <datatype> ("z_velocity"));
			datatype *pos = &(virt->index <datatype> ("z"));
			datatype value = 0.0;
			for (int j = 0; j < virt->dims ["z"] [1] - 1; ++j) {
				value += -(vel [j] - vel [j + 1]) * (vel [j] - vel [j + 1]) / ((pos [j] - pos [j + 1]) * (pos [j] - pos [j + 1])) + 1.0e-8 / ((pos [j] - pos [j + 1]) * (pos [j] - pos [j + 1]));
			}
			DEBUG (value);
			return value;
		}
	};
} /* pisces */

#endif /* end of include guard: VARDIFF_TWO_D_HPP_OQ800X4X */
