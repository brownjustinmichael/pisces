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
		int id = process_messenger.get_id ();
		int n_elements = process_messenger.get_np ();
		
		logger::log_config::set_severity (3);
		int num = 0;
		// logger::log_config::configure (&num, NULL, id, "process_%d.log");
		formats::ascii::print_headers = false;
		
		io::parameters parameters;
		
		parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
		parameters ["output.files.compare_%02i.output"] = true;
		parameters ["output.files.compare_%02i.stat"] = true;
		parameters ["output.files.compare_%02i.every"] = 10;
		parameters ["dump.file"] = "";
		parameters ["output.output"] = false;
		parameters ["time.steps"] = 100;
	
		int m = parameters.get <int> ("grid.z.points") / n_elements + 1;
		m += m % 2;

		std::vector <double> positions (n_elements + 1);
		for (int i = 0; i < n_elements + 1; ++i) {
			positions [i] = -parameters.get <double> ("grid.z.width") / 2.0 + parameters.get <double> ("grid.z.width") / n_elements * i;
		}

		int name = id;

		int n = parameters.get <int> ("grid.x.points");

		grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
		grids::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

		data::thermo_compositional_data <double> data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);
		
		std::shared_ptr <pisces::element <double>> element (new pisces::boussinesq_element <double> (horizontal_axis, vertical_axis, name, parameters, data, &process_messenger, 0x00));
		
		int n_steps = 0;
		while (n_steps < parameters.get <int> ("time.steps")) {
			element->run (n_steps, parameters.get <int> ("time.steps"), parameters.get <int> ("grid.rezone.check_every"));
		}
		
		std::string command = "./compare_cdf.py ";
		
		std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + std::string ("compare_%02i.cdf");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		command += buffer;
		
		file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("output.directory") + std::string ("compare_%02i.cdf");
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		command += " ";
		command += buffer;
		
		int exit_status = system (command.c_str ());
		
		if (exit_status != 0) {
			TS_ASSERT (false);
		}
	}
};
