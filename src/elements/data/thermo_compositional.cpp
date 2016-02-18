/*!**********************************************************************
 * \file thermo_compositional.cpp
 * /Users/justinbrown/Dropbox/pisces/src/elements/data/thermo_compositional.cpp
 ************************************************************************/

#include "thermo_compositional.hpp"

namespace data
{
	implemented_data::registrar <thermo_compositional_data> thermo_compositional_registrar ("thermo_compositional");

	std::shared_ptr <data> thermo_compositional_data::instance (grids::axis* i_axis_n, grids::axis* i_axis_m, int i_name, int n_elements, io::parameters& i_params) {
		return std::shared_ptr <data> (new thermo_compositional_data (i_axis_n, i_axis_m, i_name, n_elements, i_params));
	}

	thermo_compositional_data::thermo_compositional_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int id, int n_elements, io::parameters& i_params) : implemented_data (i_axis_n, i_axis_m, i_params, id, i_params.get <std::string> ("dump.file"), i_params.get <std::string> ("root") + i_params.get <std::string> ("dump.directory"), i_params.get <int> ("dump.every")) {
		TRACE ("Initializing...");

		initialize ("pressure", corrector);
		initialize ("composition");
		initialize ("temperature");
		initialize ("x_velocity");
		initialize ("z_velocity");
		// initialize ("density");

		// Set up the data from the input file in params
		
		if (!i_params ["output.output"].as<bool>()) {
			return;
		}
		 
		TRACE ("Setting up from input...");
		i_params ["input.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["input.directory"].as <std::string> ();
		this->setup_from <formats::netcdf> (i_params ["input"]);
		
		i_params ["output.directory"] = i_params ["root"].as <std::string> () + "/" + i_params ["output.directory"].as <std::string> ();

		DEBUG ("Setting up cartesian output...");
		std::shared_ptr <io::output> stream = this->setup_output_from <formats::netcdf> (i_params ["output.cart"]);
		// stream->append <double> ("div_u", std::shared_ptr <functors::div_functor <double>> (new functors::div_functor <double> ((*this) ["x"].ptr (), (*this) ["z"].ptr (), (*this) ["x_velocity"].ptr (), (*this) ["z_velocity"].ptr (), n, m)));

		DEBUG ("Setting up transformed output...");
		stream = this->setup_output_from <formats::netcdf> (i_params ["output.trans"], real_spectral);
		// stream->append <double> ("div_u", std::shared_ptr <functors::transform_div_functor <double>> (new functors::transform_div_functor <double> ((*this) ["x"].ptr (), (*this) ["z"].ptr (), (*this) ["x_velocity"].ptr (), (*this) ["z_velocity"].ptr (), n, m)));

		DEBUG ("Setting up profile output...");
		stream = this->setup_output_from <formats::netcdf> (i_params ["output.prof"], real_real, profile);
		for (data::iterator iter = this->begin (); iter != this->end (); ++iter) {
			std::string variable = *iter;
			DEBUG("ADDING FLUX PROF " << variable);
			DEBUG(" ");
			auto test = ((*this)[variable]);
			DEBUG(" ");
			auto test2 = ((*this)["z_velocity"]);

			stream->append <double> ("flux_" + variable, this->output_flux_prof(variable, "z_velocity"), formats::one_d);
			DEBUG("AFTER");
		}

		DEBUG ("Setting up stat output...");
		std::shared_ptr <io::output> stat_stream =
		this->setup_output_from <formats::netcdf> (i_params ["output.stat"], real_real, no_variables);

		if (!stat_stream) return;

		for (data::iterator iter = this->begin (); iter != this->end (); ++iter) {
			// For each data variable, output z_flux, average derivative across the center, average and max
			std::string variable = *iter;
			stat_stream->append <double> ("max_" + variable, this->output_max (variable), formats::scalar);
			stat_stream->append <double> ("avg_" + variable, this->output_avg (variable), formats::scalar);
			stat_stream->append <double> ("deriv_" + variable, this->output_deriv (variable), formats::scalar);
			stat_stream->append <double> ("flux_" + variable, this->output_flux (variable, "z_velocity"), formats::scalar);
			stat_stream->append <double> ("avg_flux_" + variable, this->output_avg_flux (variable, "z_velocity"), formats::scalar);
		}
	}
} /* data */
