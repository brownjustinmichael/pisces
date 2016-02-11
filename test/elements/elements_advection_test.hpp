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
#include "elements/data/thermo_compositional.hpp"
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
		
		logger::log_config::set_severity (3);
		formats::ascii::print_headers = false;
		
		io::parameters parameters(std::string (PISCES_ROOT) + "/test/elements/config.yaml");
		
		parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
		parameters ["output.cart.file"] = "advection_%02i";
		parameters ["output.cart.every"] = 10;
		parameters ["output.output"] = false;
		parameters ["dump.file"] = "";
		parameters ["time.stop"] = 1.0;
		parameters ["time.steps"] = 1000000;
		parameters ["time.max"] = 0.0001;
		parameters ["time.init"] = 0.0001;
		
		parameters ["input.file"] = "";
		
		parameters ["equations.x_velocity.ignore"] = true;
		parameters ["equations.z_velocity.ignore"] = true;
		parameters ["equations.pressure.ignore"] = true;
		parameters ["equations.composition.ignore"] = true;
		parameters ["equations.z_velocity.sources.temperature"] = 0.0;
		// parameters ["equations.z_velocity.sources.composition"] = 0.0;
		
		parameters ["equations.temperature.diffusion"] = 0.1;
		
		parameters ["grid.x.width"] = 20.0;
		parameters ["grid.z.width"] = 20.0;

		int m = 50;
		int name = id;
		int n = 300;
		double scale = 1.0 / 4.0 / acos (-1.0) / 0.1;
		double width = 2.0 * sqrt (0.1);

		grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
		grids::axis vertical_axis (m, -parameters.get <double> ("grid.z.width") / 2.0, parameters.get <double> ("grid.z.width") / 2.0, id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

		data::thermo_compositional_data data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);

		data.initialize ("temperature_0");

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data ["temperature"] [i * m + j] = data ["composition"] [i * m + j] = scale * exp (-(data ("x") [i * m + j] * data ("x") [i * m + j] + data ("z") [i * m + j] * data ("z") [i * m + j]) / width / width);
				data ["temperature_0"] [i * m + j] = scale / 2.0 * exp (-(data ("x") [i * m + j] * data ("x") [i * m + j] + data ("z") [i * m + j] * data ("z") [i * m + j]) / width / width / 2.0);
				data ["x_velocity"] [i * m + j] = 20.0;
			}
		}

		std::shared_ptr <pisces::element> element (new pisces::boussinesq_element (horizontal_axis, vertical_axis, name, parameters, data, &*process_messenger, 0x00));

		element->run ();

		double total = 0.0;
		double diff = 0.0;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				diff += (data ("temperature") [i * m + j] - data ("temperature_0") [i * m + j]) * (data ("temperature") [i * m + j] - data ("temperature_0") [i * m + j]);
				total += data ("temperature_0") [i * m + j] * data ("temperature_0") [i * m + j];	
			}
		}

		INFO ("L2 relative error is " << diff / total);
		TS_ASSERT (fabs (diff / total) < 3.0e-2);
	}
};
