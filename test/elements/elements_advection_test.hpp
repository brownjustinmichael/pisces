/*!**********************************************************************
 * \file element_test.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2015-01-28.
 * Copyright 2015 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>
#include <fstream>

#include "mpi/messenger.hpp"
#include "io/input.hpp"
#include "io/parameters.hpp"
#include "plans/plan.hpp"
#include "elements/boussinesq.hpp"
#include "logger/logger.hpp"
#include "io/formats/ascii.hpp"

class element_test_suite : public CxxTest::TestSuite
{
	std::shared_ptr <mpi::messenger> process_messenger;
		
public:
	void test_elements () {
		if (!process_messenger) process_messenger.reset (new mpi::messenger ());	

		int id = 0;
		int n_elements = 1;
		
		logger::log_config::set_severity (2);
		formats::ascii::print_headers = false;
		
		io::parameters parameters;
		
		parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
		parameters ["output.cart.file"] = "advection_%02i";
		parameters ["output.cart.every"] = 10;
		parameters ["dump.file"] = "";
		parameters ["time.stop"] = 20.0;
		parameters ["time.steps"] = 1000000;
		parameters ["time.max"] = 0.01;
		parameters ["time.init"] = 0.01;
		
		parameters ["input.file"] = "";
		
		parameters ["equations.x_velocity.ignore"] = true;
		parameters ["equations.z_velocity.ignore"] = true;
		parameters ["equations.pressure.ignore"] = true;
		parameters ["equations.composition.ignore"] = true;
		parameters ["equations.z_velocity.sources.temperature"] = 0.0;
		// parameters ["equations.z_velocity.sources.composition"] = 0.0;
		
		parameters ["equations.temperature.diffusion"] = 0.0;
		parameters ["equations.composition.diffusion"] = 0.0;
		parameters ["equations.composition.advection"] = 0.0;
		
		parameters ["grid.x.width"] = 20.0;
		parameters ["grid.z.width"] = 20.0;

		parameters ["time.mult"] = 1.02;

		int m = 200;
		int name = id;
		int n = 300;
		double scale = 1.0 / 2.0 / 3.14159;
		double width = 1.0;

		grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
		grids::axis vertical_axis (m, -parameters.get <double> ("grid.z.width") / 2.0, parameters.get <double> ("grid.z.width") / 2.0, id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

		data::thermo_compositional_data <double> data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);

		data.initialize ("temperature_0");

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data ["temperature_0"] [i * m + j] = data ["temperature"] [i * m + j] = data ["composition"] [i * m + j] = scale * exp (-(data ("x") [i * m + j] * data ("x") [i * m + j] + data ("z") [i * m + j] * data ("z") [i * m + j]) / 2.0 / width / width);
				data ["x_velocity"] [i * m + j] = 1.0;
			}
		}

		std::shared_ptr <pisces::element <double>> element (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, parameters, data, &*process_messenger, 0x00));

		element->run ();

		double total = 0.0;
		double diff = 0.0;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				diff += data ("temperature") [i * m + j] - data ("temperature_0") [i * m + j];
				total += data ("temperature_0") [i * m + j];	
			}
		}

		INFO ("L1 relative error is " << diff / total);
		TS_ASSERT (diff / total < 2.2e-4);

		// std::string command = "./compare_cdf.py ";
		
		// std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + std::string ("advection_%02i.cdf");
		// char buffer [file_format.size () * 2];
		// snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		// command += buffer;
		
		// file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("output.directory") + std::string ("advection_%02i.cdf");
		// snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		// command += " ";
		// command += buffer;
		
		// int exit_status = system (command.c_str ());
		
		// if (exit_status != 0) {
		// 	TS_ASSERT (false);
		// }
	}
};
