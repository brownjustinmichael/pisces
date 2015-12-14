/*!**********************************************************************
 * \file pseudo_element.hpp
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
	class pseudo_element : public implemented_element
	{
	protected:
		using implemented_element::element_flags;
		using implemented_element::params;
		using implemented_element::initialize;
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
		double *x_ptr; //!< The pointer to the x data, for speed
		double *z_ptr; //!< The pointer to the z data, for speed
		double *x_vel_ptr; //!< The pointer to the x velocity data, for speed
		double *z_vel_ptr; //!< The pointer to the z velocity data, for speed
		
	protected:
		using implemented_element::data;
		using implemented_element::grids;
		using implemented_element::n;
		using implemented_element::m;
		using implemented_element::equations;
		using implemented_element::messenger_ptr;
		
	public:
		using element::ptr;
		
		/*!**********************************************************************
		 * \copydoc implemented_element::implemented_element
		 ************************************************************************/
		pseudo_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~pseudo_element () {}

		static std::shared_ptr <element> instance (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) {
			return std::shared_ptr <element> (new pseudo_element (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags));
		}

		/**
		 * @copydoc element::class_name
		 */
		virtual std::string class_name() {
			return "pseudo";
		}
		
		/*!**********************************************************************
		 * \copydoc implemented_element::calculate_timestep
		 ************************************************************************/
		double calculate_timestep (int i, int j, formats::virtual_file *virtual_file = NULL);
		
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
					virtual_output->append <double> ("z", std::shared_ptr <functors::functor> (new functors::average_functor <double> (ptr ("z"), n, m)));
					virtual_output->append <double> ("z_velocity", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> (ptr ("z_velocity"), n, m)));
					virtual_output->append <double> ("temperature", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> (ptr ("temperature"), n, m)));
					virtual_output->append <double> ("composition", std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> (ptr ("composition"), n, m)));
				} else {
					FATAL ("HAVEN'T GOT A TREATMENT FOR THIS YET");
					throw 0;
					data.setup_profile (virtual_output, data::no_save);
				}
			} else {
				virtual_output.reset (new io::formatted_output <formats::virtual_format> (formats::data_grid::two_d (n, m), "two_d/pseudo/virtual_file", formats::replace_file));
				if (flags & timestep_only) {
					// If only the timestep is needed, grab the positions and velocities
					virtual_output->append <double> ("z", ptr ("z"));
					virtual_output->append <double> ("x", ptr ("x"));
					virtual_output->append <double> ("z_velocity", ptr ("z_velocity"));
					virtual_output->append <double> ("x_velocity", ptr ("x_velocity"));
					virtual_output->append <double> ("temperature", ptr ("temperature"));
					virtual_output->append <double> ("composition", ptr ("composition"));
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
		virtual formats::virtual_file *make_rezoned_virtual_file (double *positions, formats::virtual_file *virtual_file_ptr, int flags = 0x00) {
			grids::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <grids::grid> vertical_grid = implemented_element::generate_grid (&vertical_axis);
			
			pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, virtual_file_ptr, &formats::virtual_files ["two_d/pseudo/new_virtual_file"], value_buffer, inter_buffer);
			
			return &formats::virtual_files ["two_d/pseudo/new_virtual_file"];
		}
	};
} /* pisces */

#endif /* end of include guard: pseudo_TWO_D_HPP_OQ800X4X */
