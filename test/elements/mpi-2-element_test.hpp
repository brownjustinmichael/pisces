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
		io::formats::ascii::print_headers = false;
		
		io::parameters parameters;
		
		parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
		parameters ["output.stat.file"] = "compare_%02i";
		parameters ["output.stat.every"] = 10;
		parameters ["output.transform.file"] = "";
		parameters ["dump.file"] = "";
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
		
		
		std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + std::string ("compare_%02i.dat");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		std::ifstream template_stream (buffer, std::ifstream::in);
		
		file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("output.directory") + std::string ("compare_%02i.dat");
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);
		std::ifstream compare_stream (buffer, std::ifstream::in);
		
		double template_string, compare_string;
		
		while (!template_stream.eof ()) {
			template_stream >> template_string;
			compare_stream >> compare_string;
			
			TS_ASSERT_DELTA (template_string, compare_string, 1.0E-8);
		}
	}
};
