/*!**********************************************************************
 * \file pseudo.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef pseudo_TWO_D_HPP_OQ800X4X
#define pseudo_TWO_D_HPP_OQ800X4X

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
	template <class datatype>
	class pseudo_element : public implemented_element <datatype>
	{
	protected:
		using implemented_element <datatype>::element_flags;
		using implemented_element <datatype>::params;
		using implemented_element <datatype>::initialize;
		using implemented_element <datatype>::name;
		using implemented_element <datatype>::cell_n;
		using implemented_element <datatype>::cell_m;
		using implemented_element <datatype>::matrix_ptr;
		using implemented_element <datatype>::timestep;
		using implemented_element <datatype>::duration;
		using implemented_element <datatype>::value_buffer;
		using implemented_element <datatype>::inter_buffer;
		
		std::vector <datatype> diffusion; //!< A vector of diffusion data, for background diffusion
		datatype advection_coeff; //!< The advection coefficient, for speed
		datatype cfl; //!< The cfl coefficient, for speed
		datatype *x_ptr; //!< The pointer to the x data, for speed
		datatype *z_ptr; //!< The pointer to the z data, for speed
		datatype *x_vel_ptr; //!< The pointer to the x velocity data, for speed
		datatype *z_vel_ptr; //!< The pointer to the z velocity data, for speed
		
	protected:
		using implemented_element <datatype>::data;
		using implemented_element <datatype>::grids;
		using implemented_element <datatype>::n;
		using implemented_element <datatype>::m;
		using implemented_element <datatype>::equations;
		using implemented_element <datatype>::messenger_ptr;
		
	public:
		using element <datatype>::ptr;
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		pseudo_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~pseudo_element () {}
		
		/*!**********************************************************************
		 * \copydoc implemented_element::calculate_timestep
		 ************************************************************************/
		datatype calculate_timestep (int i, int j, formats::virtual_file *virtual_file = NULL);
		
		/*!**********************************************************************
		 * \copydoc implemented_element::make_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_virtual_file (int flags = 0x00) {
			std::shared_ptr <io::output> virtual_output;
			if (flags & profile_only) {
				// If only the profile is desired, just build that
				virtual_output.reset (new io::formatted_output <formats::virtual_format> (formats::data_grid::two_d (1, m), "two_d/pseudo/virtual_file", formats::replace_file));
				if (flags & timestep_only) {
					// If only the timestep is needed, only load z and x_velocity
					virtual_output->append <datatype> ("z", std::shared_ptr <functors::functor> (new functors::average_functor <datatype> (ptr ("z"), n, m)));
					virtual_output->append <datatype> ("z_velocity", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <datatype> (ptr ("z_velocity"), n, m)));
					virtual_output->append <datatype> ("temperature", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <datatype> (ptr ("temperature"), n, m)));
					virtual_output->append <datatype> ("composition", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <datatype> (ptr ("composition"), n, m)));
				} else {
					FATAL ("HAVEN'T GOT A TREATMENT FOR THIS YET");
					throw 0;
					data.setup_profile (virtual_output, data::no_save);
				}
			} else {
				virtual_output.reset (new io::formatted_output <formats::virtual_format> (formats::data_grid::two_d (n, m), "two_d/pseudo/virtual_file", formats::replace_file));
				if (flags & timestep_only) {
					// If only the timestep is needed, grab the positions and velocities
					virtual_output->append <datatype> ("z", ptr ("z"));
					virtual_output->append <datatype> ("x", ptr ("x"));
					virtual_output->append <datatype> ("z_velocity", ptr ("z_velocity"));
					virtual_output->append <datatype> ("x_velocity", ptr ("x_velocity"));
					virtual_output->append <datatype> ("temperature", ptr ("temperature"));
					virtual_output->append <datatype> ("composition", ptr ("composition"));
				} else {
					// Load the whole dataset
					data.setup_output (virtual_output, data::no_save);
				}
			}
			virtual_output->to_file ();
			return &formats::virtual_files ["two_d/pseudo/virtual_file"];
		}
		
		/*!**********************************************************************
		 * \copydoc implemented_element::make_rezoned_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_rezoned_virtual_file (datatype *positions, formats::virtual_file *virtual_file_ptr, int flags = 0x00) {
			grids::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <grids::grid <datatype>> vertical_grid = implemented_element <datatype>::generate_grid (&vertical_axis);
			
			pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, virtual_file_ptr, &formats::virtual_files ["two_d/pseudo/new_virtual_file"], value_buffer, inter_buffer);
			
			return &formats::virtual_files ["two_d/pseudo/new_virtual_file"];
		}
	};
} /* pisces */

#endif /* end of include guard: pseudo_TWO_D_HPP_OQ800X4X */
