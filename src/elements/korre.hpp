/*!**********************************************************************
 * \file korre.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef KORRE_TWO_D_HPP_OQ800X4X
#define KORRE_TWO_D_HPP_OQ800X4X

#include "mpi/messenger.hpp"

#include "io/functors/div.hpp"
#include "io/functors/product.hpp"

#include "boussinesq.hpp"

namespace pisces
{
	/*!**********************************************************************
	 * \brief An element built to use the Boussinesq approximation
	 ************************************************************************/
	class korre_element : public boussinesq_element
	{
	protected:
		using boussinesq_element::element_flags;
		using boussinesq_element::params;
		using boussinesq_element::name;
		using boussinesq_element::cell_n;
		using boussinesq_element::cell_m;
		using boussinesq_element::matrix_ptr;
		using boussinesq_element::timestep;
		using boussinesq_element::duration;
		using boussinesq_element::value_buffer;
		using boussinesq_element::inter_buffer;
		
	protected:
		using boussinesq_element::data;
		using boussinesq_element::grids;
		using boussinesq_element::n;
		using boussinesq_element::m;
		using boussinesq_element::equations;
		using boussinesq_element::messenger_ptr;
		
	public:
		using boussinesq_element::initialize;

		using element::ptr;
		
		/*!**********************************************************************
		 * \copydoc boussinesq_element::boussinesq_element
		 ************************************************************************/
		korre_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~korre_element () {}

		static std::shared_ptr <element> instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) {
			return std::shared_ptr <element> (new korre_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags));
		}

		virtual std::string class_name() {
			return "korre";
		}
	};
} /* pisces */

#endif /* end of include guard: KORRE_TWO_D_HPP_OQ800X4X */
