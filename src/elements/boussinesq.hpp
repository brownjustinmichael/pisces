/*!**********************************************************************
 * \file boussinesq.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUSSINESQ_TWO_D_HPP_OQ800X4X
#define BOUSSINESQ_TWO_D_HPP_OQ800X4X

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
	class boussinesq_element : public implemented_element
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
		using implemented_element::value_buffer;
		using implemented_element::inter_buffer;
		
		std::vector <double> diffusion; //!< A vector of diffusion data, for background diffusion
		double advection_coeff; //!< The advection coefficient, for speed
		double cfl; //!< The cfl coefficient, for speed
		double allow; //!< The change allowance, for speed
		double *x_ptr; //!< The pointer to the x data, for speed
		double *z_ptr; //!< The pointer to the z data, for speed
		double *x_vel_ptr; //!< The pointer to the x velocity data, for speed
		double *z_vel_ptr; //!< The pointer to the z velocity data, for speed
		double *t_spectral; //!< A pointer to the T spectral data, for speed
		double *t_rhs; //!< A pointer to the T rhs data, for speed
		
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
		boussinesq_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~boussinesq_element () {}

		static std::shared_ptr <element> instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);

		virtual std::string class_name();
		
		/*!**********************************************************************
		 * \copydoc implemented_element::calculate_timestep
		 ************************************************************************/
		double calculate_timestep (int i, int j, formats::virtual_file *virtual_file = NULL);
		
		/*!**********************************************************************
		 * \copydoc implemented_element::make_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_virtual_file (int flags = 0x00);
		
		/*!**********************************************************************
		 * \copydoc implemented_element::make_rezoned_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_rezoned_virtual_file (double *positions, formats::virtual_file *virtual_file_ptr, int flags = 0x00);
	};
} /* pisces */

namespace data
{
	/*!**********************************************************************
	 * \brief A data object designed to hold and output thermo-compositional data
	 ************************************************************************/
	class thermo_compositional_data : public implemented_data
	{
	protected:
		using implemented_data::n;
		using implemented_data::m;
		using implemented_data::grid_m;
		using implemented_data::iterator;
		using implemented_data::duration;
		
		std::vector <double> area; //!< A vector containing the area of each cell, for weighted averages
		
	public:
		using implemented_data::initialize;

		/*!**********************************************************************
		 * \param i_axis_n The horizontal axis object
		 * \param i_axis_m The vertical axis object
		 * \param i_name The integer name of the element
		 * \param n_elements The total number of elements
		 * \param i_params The parameters object associated with the run
		 ************************************************************************/
		thermo_compositional_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int i_name, int n_elements, io::parameters& i_params);
		
		virtual ~thermo_compositional_data () {}
	};
} /* data */

#endif /* end of include guard: BOUSSINESQ_TWO_D_HPP_OQ800X4X */
