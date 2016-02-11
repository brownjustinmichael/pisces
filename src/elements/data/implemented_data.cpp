/*!**********************************************************************
 * \file implemented_data.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-09.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/


#include "implemented_data.hpp"

#include "io/functors/profile.hpp"

namespace data
{
	std::map <std::string, implemented_data::data_function> & implemented_data::registry()
	{
	    static std::map <std::string, implemented_data::data_function> impl;
	    return impl;
	}

	implemented_data::implemented_data (grids::axis *i_axis_n, grids::axis *i_axis_m, io::parameters &i_params, int i_name, std::string dump_file, std::string dump_directory, int dump_every) : 
	data (i_params, i_name), 
	grid_n (std::shared_ptr <grids::grid> (new grids::horizontal::grid (i_axis_n))), 
	grid_m (std::shared_ptr <grids::grid> (new grids::vertical::grid (i_axis_m))), 
	n (grid_n->get_n ()), 
	m (grid_m->get_n ()) {
		// Set up output
		const formats::data_grid o_grid = formats::data_grid::two_d (n, m, 0, 0, 0, 0);
		
		std::shared_ptr <io::output> dump_stream;
		if (dump_file != "") {
			std::string file_format = dump_directory + dump_file;
			char buffer [file_format.size () * 2];
			snprintf (buffer, file_format.size () * 2, file_format.c_str (), i_name);

			dump_stream.reset (new io::replace_output <formats::netcdf> (o_grid, buffer, dump_every));
			this->setup_dump (dump_stream);
		}
		
		this->initialize ("x");
		this->initialize ("z");

		// For weighted averages, calculate area
		area.resize (n * m);
		for (int i = 1; i < n; ++i) {
			for (int j = 1; j < m; ++j) {
				area [i * m + j] = ((*(this->grid_n)) [i] - (*(this->grid_n)) [i - 1]) * ((*(this->grid_m)) [j] - (*(this->grid_m)) [j - 1]);
			}
		}

		weights = &area [0];
		
		/*
			TODO Clean dump file generation
		*/
	}

	std::shared_ptr <data> implemented_data::instance (std::string const & name, grids::axis* i_axis_n, grids::axis* i_axis_m, int i_name, int n_elements, io::parameters& i_params) {
		auto it = registry ().find (name);
		data_function fp = (it == registry ().end () ? NULL : (it->second));
		if (fp) {
			return fp (i_axis_n, i_axis_m, i_name, n_elements, i_params);
		} else {
			WARN("No data type associated with name " << name);
			return std::shared_ptr <data> ();
		}
	}

	void implemented_data::setup_profile (std::shared_ptr <io::output> output_ptr, int flags) {
		typedef std::map <std::string, std::shared_ptr <grids::variable>>::iterator iterator;
		for (iterator iter = data::variables.begin (); iter != data::variables.end (); ++iter) {
			output_ptr->append <double> (iter->first, (*this) (iter->first));
			output_ptr->append <double> ("rms_" + iter->first, std::shared_ptr <functors::functor> (new functors::root_mean_square_functor <double> ((*this) (iter->first), n, m)));
			output_ptr->append <double> ("avg_" + iter->first, std::shared_ptr <functors::functor> (new functors::average_functor <double> ((*this) (iter->first), n, m)));
		}
		output_ptr->append <double> ("t", &(duration), formats::scalar);
		output_ptr->append <const int> ("mode", &(data::get_mode ()), formats::scalar);

		data::setup_profile (output_ptr, flags);
	}

	std::shared_ptr <functors::functor> implemented_data::output_max (std::string variable) {
		return std::shared_ptr <functors::functor> (new functors::max_functor <double> (n, m, (*this) (variable)));
	}

	std::shared_ptr <functors::functor> implemented_data::output_avg (std::string variable) {
		return std::shared_ptr <functors::functor> (new functors::weighted_average_functor <double> (n, m, this->weights, (*this) (variable)));
	}

	std::shared_ptr <functors::functor> implemented_data::output_deriv (std::string variable, int slice_index) {
		return std::shared_ptr <functors::functor> (new functors::average_functor <double> (n, 1, std::shared_ptr <functors::functor> (new functors::slice_functor <double> (n, m, slice_index < 0 ? m / 2 : slice_index, std::shared_ptr <functors::functor> (new functors::deriv_functor <double> ((*this) (variable), n, m, &(*grid_m) [0]))))));
	}

	std::shared_ptr <functors::functor> implemented_data::output_flux (std::string variable, std::string velocity, int slice_index) {
		return std::shared_ptr <functors::functor> (new functors::average_functor <double> (n, 1, std::shared_ptr <functors::functor> (new functors::slice_functor <double> (n, m, slice_index < 0 ? m / 2 : slice_index, std::shared_ptr <functors::functor> (new functors::product_functor <double> (n, m, (*this) (velocity), (*this) (variable)))))));
	}

	std::shared_ptr <functors::functor> implemented_data::output_avg_flux (std::string variable, std::string velocity) {
		return std::shared_ptr <functors::functor> (new functors::weighted_average_functor <double> (n, m, this->weights, std::shared_ptr <functors::functor> (new functors::product_functor <double> (n, m, (*this) (velocity), (*this) (variable)))));
	}

	std::shared_ptr <functors::functor> implemented_data::output_prof (std::string variable) {
		return std::shared_ptr <functors::functor> (new functors::profile_functor <double> (n, m, (*this) (variable)));
	}

	std::shared_ptr <functors::functor> implemented_data::output_deriv_prof (std::string variable) {
		return std::shared_ptr <functors::functor> (new functors::profile_functor <double> (n, m, std::shared_ptr <functors::functor> (new functors::deriv_functor <double> ((*this) (variable), n, m, &(*grid_m) [0]))));
	}

	double *implemented_data::operator() (const std::string &name, int state, int i, int j) {
		return data::operator() (name, state, i * m + j);
	}

	grids::variable &implemented_data::_initialize (std::string name, int i_flags) {
		TRACE ("Initializing " << name << "...");
		// Size allowing for real FFT buffer
		if (i_flags & uniform_n) {
			data::variables [name] = std::shared_ptr <grids::variable> (new grids::variable (*grid_m, flags, name));
			return *variables [name];
		} else if (i_flags & uniform_m) {
			data::variables [name] = std::shared_ptr <grids::variable> (new grids::variable (*grid_n, flags, name));
			return *variables [name];
		}
		data::variables [name] = std::shared_ptr <grids::variable> (new grids::variable (*grid_n, *grid_m, flags, name, 3, i_flags & vector2D ? 2 : (i_flags & vector3D ? 3 : 1)));
		if (name == "x") {
			for (int j = 0; j < m; ++j) {
				linalg::copy (n, &((*grid_n) [0]), (*this) (name, real_real, 0, j), 1, m);
			}
			transformers [name] = std::shared_ptr <plans::transforms::transformer> ();
		} else if (name == "z") {
			for (int i = 0; i < n; ++i) {
				linalg::copy (m, &((*grid_m) [0]), (*this) (name, real_real, i));
			}
			transformers [name] = std::shared_ptr <plans::transforms::transformer> ();
		} else {
			transformers [name] = std::shared_ptr <plans::transforms::transformer > (new plans::transforms::implemented_transformer ((*this) [name], &(flags), &((*this) [name].component_flags)));
		}
		TRACE ("Done.");
		return *variables [name];
	}
}
