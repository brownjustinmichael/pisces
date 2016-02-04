/*!**********************************************************************
 * \file diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_TWO_D_HPP_OQ800X4X
#define DIFFUSION_TWO_D_HPP_OQ800X4X

#include "mpi/messenger.hpp"

#include "io/functors/div.hpp"
#include "io/functors/product.hpp"

#include "implemented_element.hpp"

#include "data/data.hpp"

namespace pisces
{
	/*!**********************************************************************
	 * \brief An element built to use the Boussinesq approximation
	 ************************************************************************/
	class diffusion_element : public implemented_element
	{
	protected:
		using implemented_element::element_flags;
		using implemented_element::params;
		using implemented_element::name;
		using implemented_element::cell_n;
		using implemented_element::cell_m;
		using implemented_element::matrix_ptr;
		using implemented_element::timestep;
		using implemented_element::duration;
		
	protected:
		using implemented_element::data;
		using implemented_element::grids;
		using implemented_element::n;
		using implemented_element::m;
		using implemented_element::equations;
		using implemented_element::messenger_ptr;
		
	public:
		using implemented_element::initialize;
		using element::ptr;
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		diffusion_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~diffusion_element () {}

		static std::shared_ptr <element> instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);

		virtual std::string class_name();
		
		/*!**********************************************************************
		 * \copydoc implemented_element::calculate_timestep
		 ************************************************************************/
		double calculate_timestep (int i, int j, formats::virtual_file *virtual_file = NULL);
	};
} /* pisces */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_OQ800X4X */
