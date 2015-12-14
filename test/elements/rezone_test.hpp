/*!**********************************************************************
 * \file rezone_test.hpp
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

// double rezone_merit (pisces::element *, formats::virtual_file *) {
	
// }


class rezone_test_suite : public CxxTest::TestSuite
{
	std::shared_ptr <mpi::messenger> process_messenger;
	
public:
	void test_elements () {
	// 	if (!process_messenger) process_messenger.reset (new mpi::messenger ());
	// 	int id = 0;
	// 	int n_elements = 1;
		
	// 	logger::log_config::set_severity (3);
	// 	formats::ascii::print_headers = false;
		
	// 	io::parameters parameters;
		
	// 	parameters ["root"] = std::string (PISCES_ROOT) + "/test/elements/";
	// 	parameters ["output.cart.file"] = "rezone_%02i";
	// 	parameters ["output.cart.every"] = 100;
	// 	parameters ["output.trans.file"] = "rezone_t_%02i";
	// 	parameters ["output.trans.every"] = 100;

	// 	// parameters ["output.output"] = false;
	// 	parameters ["dump.file"] = "";
		
	// 	parameters ["time.stop"] = 1.0;
	// 	parameters ["time.max"] = 1.e-2;
	// 	parameters ["time.init"] = 1.e-6;
	// 	parameters ["time.mult"] = 1.01;
	// 	parameters ["time.steps"] = 10000;
		
	// 	parameters ["input.file"] = "";
		
	// 	parameters ["equations.composition.ignore"] = true;
		
	// 	parameters ["equations.temperature.diffusion"] = 0.0;
	// 	// parameters ["equations.temperature.sources.z_velocity"] = 1.0;

	// 	parameters ["equations.temperature.top.value"] = 0.0;
	// 	parameters ["equations.temperature.bottom.value"] = 1.0;

	// 	parameters ["equations.velocity.diffusion"] = 1.0;
	// 	parameters ["equations.z_velocity.sources.temperature"] = 0.0;
		
	// 	parameters ["grid.x.width"] = 2.8284;
	// 	parameters ["grid.z.width"] = 1.0;
	
	// 	int m = 100;
	// 	int name = id;
	// 	int n = 100;
	// 	double scale = 0.00001;

	// 	grids::axis horizontal_axis (n, -parameters.get <double> ("grid.x.width") / 2.0, parameters.get <double> ("grid.x.width") / 2.0);
	// 	grids::axis vertical_axis (m, -parameters.get <double> ("grid.z.width") / 2.0, parameters.get <double> ("grid.z.width") / 2.0, id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);

	// 	data::thermo_compositional_data data (&horizontal_axis, &vertical_axis, id, n_elements, parameters);

	// 	std::shared_ptr <pisces::element> element (new pisces::boussinesq_element (horizontal_axis, vertical_axis, name, parameters, data, &*process_messenger, 0x00));

	// 	for (int q = 0; q < 100; ++q)
	// 	{
	// 		for (int i = 0; i < n; ++i) {
	// 			for (int j = 0; j < m; ++j) {
	// 				if (data ["z"] [j] < q / 100.) {
	// 					data ["temperature"] = 0.0;
	// 				} else {
	// 					data ["temperature"] = 1.0;
	// 				}
	// 			}
	// 		}

	// 		data ["temperature"].update ();

	// 		formats::virtual_file *virt = element->rezone_minimize_ts (&positions [0], parameters.get <double> ("grid.rezone.min_size"), parameters.get <double> ("grid.rezone.max_size"), parameters.get <int> ("grid.rezone.n_tries"), parameters.get <int> ("grid.rezone.iters_fixed_t"), parameters.get <double> ("grid.rezone.step_size"), parameters.get <double> ("grid.rezone.k"), parameters.get <double> ("grid.rezone.t_initial"), parameters.get <double> ("grid.rezone.mu_t"), parameters.get <double> ("grid.rezone.t_min"), );
			
	// 		if (virt) {
	// 			formats::virtual_files ["main/virtual_file"] = *virt;
	// 			grids::axis vertical_axis (m, positions [id], positions [id + 1], id == 0 ? 0 : 1, id == n_elements - 1 ? 0 : 1);
			
	// 			virtual_input.reset (new io::formatted_input <formats::virtual_format> (formats::data_grid::two_d (n, m), "main/virtual_file"));
	// 			data.setup (virtual_input);
			
	// 			element = pisces::implemented_element::instance (parameters ["element"].as <std::string> (), horizontal_axis, vertical_axis, name, parameters, data, &process_messenger, 0x00);
	// 		}
	// 	}

	// 	element->run ();

	// 	TS_ASSERT (fabs ((*(double *) (avg_flux->calculate ())) - *(double *) (avg_deriv->calculate ()) - 3.04) < 0.02);
	}
};
