/*!**********************************************************************
 * \file pseudo.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef pseudo_TWO_D_HPP_OQ800X4X
#define pseudo_TWO_D_HPP_OQ800X4X

#include "boussinesq.hpp"

namespace pisces
{
	/*!**********************************************************************
	 * \brief An element built to use the pseudo approximation
	 ************************************************************************/
	template <class datatype>
	class pseudo_element : public boussinesq_element <datatype>
	{	
	protected:
		using boussinesq_element <datatype>::equations;
		using boussinesq_element <datatype>::data;
		using boussinesq_element <datatype>::m;
		using boussinesq_element <datatype>::grids;
		std::vector <datatype> pressure;

	public:
		using element <datatype>::ptr;
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		pseudo_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags, bool load_diffusion = true);
		
		virtual ~pseudo_element () {}
	};
} /* pisces */

#endif /* end of include guard: pseudo_TWO_D_HPP_OQ800X4X */
