/*!**********************************************************************
 * \file boussinesq_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-03.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUSSINESQ_TWO_D_HPP_OQ800X4X
#define BOUSSINESQ_TWO_D_HPP_OQ800X4X

#include "io/functors/div.hpp"
#include "io/functors/product.hpp"

#include "implemented_element.hpp"

namespace pisces
{
	template <class datatype>
	class boussinesq_element : public implemented_element <datatype>
	{
	public:
		boussinesq_element (plans::axis i_axis_n, plans::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags);
		
		virtual ~boussinesq_element () {}
		
		datatype calculate_timestep (int i, int j, io::formats::virtual_file *virtual_file = NULL);
		
		virtual void setup_output (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			element <datatype>::setup_output (output_ptr, flags);
			output_ptr->template append <double> ("div", new io::functors::div_functor <datatype> (ptr (x_position), ptr (z_position), ptr (x_velocity), ptr (z_velocity), n, m));
		}
		
		virtual void setup_stat (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			element <datatype>::setup_stat (output_ptr, flags);
			area.resize (n * m);
			for (int i = 1; i < n; ++i) {
				for (int j = 1; j < m; ++j) {
					area [i * m + j] = ((*(grids [0])) [i] - (*(grids [0])) [i - 1]) * ((*(grids [1])) [j] - (*(grids [1])) [j - 1]);
				}
			}
			output_ptr->template append <double> ("wT", new io::functors::weighted_average_functor <datatype> (n, m, &area [0], new io::functors::product_functor <datatype> (n, m, ptr (z_velocity), ptr (temperature))), io::scalar);
			output_ptr->template append <double> ("wS", new io::functors::weighted_average_functor <datatype> (n, m, &area [0], new io::functors::product_functor <datatype> (n, m, ptr (z_velocity), ptr (composition))), io::scalar);
		}
		
		virtual io::formats::virtual_file *make_virtual_file (int flags = 0x00) {
			std::shared_ptr <io::output> virtual_output;
			if (flags & profile_only) {
				virtual_output.reset (new io::formatted_output <io::formats::two_d::virtual_format> (io::data_grid::two_d (1, m), "two_d/boussinesq/virtual_file", io::replace_file));
				if (flags & timestep_only) {
					virtual_output->append <datatype> ("z", new io::functors::average_functor <datatype> (ptr (z_position), n, m));
					virtual_output->append <datatype> ("w", new io::functors::root_mean_square_functor <datatype> (ptr (z_velocity), n, m));
				} else {
					element <datatype>::setup_profile (virtual_output);
				}
			} else {
				virtual_output.reset (new io::formatted_output <io::formats::two_d::virtual_format> (io::data_grid::two_d (n, m), "two_d/boussinesq/virtual_file", io::replace_file));
				if (flags & timestep_only) {
					virtual_output->append <datatype> ("z", ptr (z_position));
					virtual_output->append <datatype> ("x", ptr (x_position));
					virtual_output->append <datatype> ("w", ptr (z_velocity));
					virtual_output->append <datatype> ("u", ptr (x_velocity));
				} else {
					element <datatype>::setup_output (virtual_output);
				}
			}
			virtual_output->to_file ();
			return &io::virtual_files ["two_d/boussinesq/virtual_file"];
		}
		
		virtual io::formats::virtual_file *make_rezoned_virtual_file (datatype *positions, io::formats::virtual_file *old_virtual_file, int flags = 0x00) {
			plans::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <plans::grid <datatype>> vertical_grid = implemented_element <datatype>::generate_grid (&vertical_axis);
			
			pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, old_virtual_file, &io::virtual_files ["two_d/boussinesq/new_virtual_file"]);
			
			return &io::virtual_files ["two_d/boussinesq/new_virtual_file"];
		}
		
		using element <datatype>::ptr;
	private:
		using implemented_element <datatype>::element_flags;
		using implemented_element <datatype>::params;
		using implemented_element <datatype>::initialize;
		using implemented_element <datatype>::n;
		using implemented_element <datatype>::m;
		using implemented_element <datatype>::name;
		using implemented_element <datatype>::cell_n;
		using implemented_element <datatype>::cell_m;
		using implemented_element <datatype>::grids;
		using implemented_element <datatype>::matrix_ptr;
		using implemented_element <datatype>::messenger_ptr;
		using implemented_element <datatype>::timestep;
		using implemented_element <datatype>::duration;
		using implemented_element <datatype>::alpha_0;
		using implemented_element <datatype>::alpha_n;
		using implemented_element <datatype>::solvers;
		using implemented_element <datatype>::scalars;
		using implemented_element <datatype>::scalar_names;
		using implemented_element <datatype>::data;
		
		std::vector <datatype> area;
		datatype advection_coeff, cfl, *x_ptr, *z_ptr, *x_vel_ptr, *z_vel_ptr;
	};
} /* pisces */

#endif /* end of include guard: BOUSSINESQ_TWO_D_HPP_OQ800X4X */
