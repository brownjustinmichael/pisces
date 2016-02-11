/*!**********************************************************************
 * \file scalar.cpp
 * /Users/justinbrown/Dropbox/pisces/src/elements/data/scalar.cpp
 ************************************************************************/

#include "scalar.hpp"

namespace data
{
	implemented_data::registrar <scalar_data> scalar_registrar ("scalar");

	std::shared_ptr <data> scalar_data::instance (grids::axis* i_axis_n, grids::axis* i_axis_m, int i_name, int n_elements, io::parameters& i_params) {
		return std::shared_ptr <data> (new scalar_data (i_axis_n, i_axis_m, i_name, n_elements, i_params));
	}

	scalar_data::scalar_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data (i_axis_n, i_axis_m, i_params, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
		TRACE ("Initializing...");

		initialize ("scalar");
		initialize ("x_velocity");
		initialize ("z_velocity");

		// Set up the data from the input file in params
		TRACE ("Setting up from input...");
		i_params ["input.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["input.directory"].as <std::string> ();
		this->setup_from <formats::netcdf> (i_params ["input"]);
		
		i_params ["output.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["output.directory"].as <std::string> ();

		TRACE ("Setting up cartesian output...");
		std::shared_ptr <io::output> stream = this->setup_output_from <formats::netcdf> (i_params ["output.cart"]);
	}
} /* data */
