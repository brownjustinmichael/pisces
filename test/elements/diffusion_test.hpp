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
	mpi::messenger process_messenger;
	
public:
	void test_elements () {
		int id = 0;
		int n_elements = 1;
		
		logger::log_config::set_severity (3);
		formats::ascii::print_headers = false;
		
		io::parameters parameters;
		
		parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
		parameters ["output.files.diffusion_%02i.output"] = true;
		parameters ["output.files.diffusion_%02i.every"] = 10;
		parameters ["output.files.diffusion_%02i.stat"] = true;
		parameters ["output.output"] = false;
		parameters ["dump.file"] = "";
		
		parameters ["time.steps"] = 200;
		parameters ["time.max"] = 0.1;
		parameters ["time.init"] = 0.1;
		
		parameters ["input.file"] = "";
		
		parameters ["equations.x_velocity.ignore"] = true;
		parameters ["equations.z_velocity.ignore"] = true;
		parameters ["equations.composition.ignore"] = true;
		
		parameters ["equations.temperature.advection"] = 0.0;
		
		parameters ["grid.x.width"] = 20.0;
		parameters ["grid.z.width"] = 20.0;
	
		int m = 200;
		int name = id;
		int n = 300;
		double scale = 1.0 / 2.0 / 3.14159;
		double width = 1.0;

		grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
		grids::axis vertical_axis (m, -parameters.get <double> ("grid.z.width") / 2.0, parameters.get <double> ("grid.z.width") / 2.0, id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

		data::thermo_compositional_data <double> data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data ("temperature") [i * m + j] = scale * exp (-(data ("x") [i * m + j] * data ("x") [i * m + j] + data ("z") [i * m + j] * data ("z") [i * m + j]) / 2.0 / width / width);
			}
		}

		std::shared_ptr <pisces::element <double>> element (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, parameters, data, &process_messenger, 0x00));

		int n_steps = 0;
		while (n_steps < parameters.get <int> ("time.steps")) {
			element->run (n_steps, parameters.get <int> ("time.steps"), parameters.get <int> ("grid.rezone.check_every"));
		}
		
		std::string command = "./compare_cdf.py ";
		
		std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + std::string ("diffusion_%02i.cdf");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		command += buffer;
		
		file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("output.directory") + std::string ("diffusion_%02i.cdf");
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		command += " ";
		command += buffer;
		
		int exit_status = system (command.c_str ());
		
		if (exit_status != 0) {
			TS_ASSERT (false);
		}
	}
};
