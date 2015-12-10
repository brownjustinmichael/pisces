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

class convection_test_suite : public CxxTest::TestSuite
{
	std::shared_ptr <mpi::messenger> process_messenger;
	
public:
	void test_elements () {
		if (!process_messenger) process_messenger.reset (new mpi::messenger ());
		int id = 0;
		int n_elements = 1;
		
		logger::log_config::set_severity (3);
		formats::ascii::print_headers = false;
		
		io::parameters parameters;
		
		parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
		parameters ["output.cart.file"] = "convect_%02i";
		parameters ["output.cart.every"] = 1;
		parameters ["output.trans.file"] = "convect_t_%02i";
		parameters ["output.trans.every"] = 100;

		parameters ["output.output"] = false;
		parameters ["dump.file"] = "";
		
		parameters ["time.stop"] = 1.0;
		parameters ["time.max"] = 1.e-2;
		parameters ["time.init"] = 1.e-6;
		parameters ["time.mult"] = 1.01;
		parameters ["time.steps"] = 10000;
		parameters ["time.cfl"] = 1.e6;
		
		parameters ["input.file"] = "";
		
		parameters ["equations.composition.ignore"] = true;
		
		parameters ["equations.temperature.diffusion"] = 1.0;
		parameters ["equations.temperature.advection"] = 0.0;
		parameters ["equations.temperature.sources.z_velocity"] = 1.0;

		parameters ["equations.temperature.top.value"] = 0.0;
		parameters ["equations.temperature.bottom.value"] = 0.0;

		parameters ["equations.velocity.diffusion"] = 6.8;
		parameters ["equations.velocity.advection"] = 0.0;
		parameters ["equations.z_velocity.sources.temperature"] = 4. * 657.5 * 6.8;
		parameters ["equations.z_velocity.sources.composition"] = 0.;
		
		parameters ["grid.x.width"] = 2.8285;
		parameters ["grid.z.width"] = 1.0;
	
		int m = 100;
		int name = id;
		int n = 4;
		double scale = 0.00001;

		grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
		grids::axis vertical_axis (m, -parameters.get <double> ("grid.z.width") / 2.0, parameters.get <double> ("grid.z.width") / 2.0, id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

		data::thermo_compositional_data data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);
		data.initialize ("temperature_0");

		auto avg_flux = data.output_flux ("temperature", "z_velocity");

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				data ["temperature_0"] [i * m + j] = data ["temperature"] [i * m + j] = scale * (rand () % 100 / 100.) * cos (data ["z"] [j] / acos (-1.));// + (-data ["z"] [j] + 0.5);
			}
		}

		std::shared_ptr <pisces::element> element (new pisces::boussinesq_element (horizontal_axis, vertical_axis, name, parameters, data, &*process_messenger, 0x00));

		element->run ();

		parameters ["time.stop"] = 1.1;

		double temp = data ["temperature"].ptr (real_spectral) [2 * m + 50];

		element->run ();

		INFO ("At 1: " << temp);
		INFO ("At 1.1: " << data ["temperature"].ptr (real_spectral) [2 * m + 50]);
		
		TS_ASSERT (fabs (data ["temperature"].ptr (real_spectral) [2 * m + 50] / temp - 21.3408) < 3.0);
	}
};
