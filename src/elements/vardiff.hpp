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
		
	public:
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		vardiff_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~vardiff_element () {}
	};
} /* pisces */

#endif /* end of include guard: VARDIFF_TWO_D_HPP_OQ800X4X */
