/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include "two_d/plan_two_d.hpp"
#include "bases/grid.hpp"
#include "bases/messenger.hpp"
#include "utils/io.hpp"
#include <sstream>
#include <iomanip>

using namespace two_d::fourier::cosine;

int main (int argc, char *argv[])
{
	try {
		int id = 0, n_elements = 1;

		// Initialize messenger
		bases::messenger process_messenger (&argc, &argv, 2);

		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();

		log_config::configure (&argc, &argv);
	
		io::parameters config ("../input/config.yaml");
		if (!config ["grid.x.points"].IsDefined ()) config ["grid.x.points"] = 64;
		if (!config ["grid.z.points"].IsDefined ()) config ["grid.z.points"] = 64;

		int m = config.get <int> ("grid.z.points") / n_elements + 1;
		int n = config.get <int> ("grid.x.points");
		double position_m0 = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * id;
		double position_mm = -config.get <double> ("grid.z.width") / 2.0 + config.get <double> ("grid.z.width") / n_elements * (id + 1);
		double position_n0 = -config.get <double> ("grid.x.width") / 2.0;
		double position_nn = config.get <double> ("grid.x.width") / 2.0;
	
		int excess_0;
		int excess_n;
		if (id == 0) {
			excess_0 = 0;
		} else {
			excess_0 = 1;
		}
		if (id == n_elements - 1) {
			excess_n = 0;
		} else {
			excess_n = 1;
		}
	
		horizontal_grid <double>::type horizontal_grid (new bases::axis (n, position_n0, position_nn));
		vertical_grid <double>::type vertical_grid (new bases::axis (m, position_m0, position_mm, excess_0, excess_n));

		std::vector <double> temp (n * m);

		double scale = config.get <double> ("init.scale");
		double mean = config.get <double> ("init.mean");
		double sigma = config.get <double> ("init.sigma");
		double width = config.get <double> ("grid.z.width");
		double height;
		srand (id);
		height = std::max (scale * std::exp (- (-width / 2.0 - mean) * (-width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma));
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				temp [i * m + j] = scale * (std::exp (- (horizontal_grid [i] - mean) * (horizontal_grid [i] - mean) / 2.0 / sigma / sigma) - height) * (std::exp (- (vertical_grid [j] - mean) * (vertical_grid [j] - mean) / 2.0 / sigma / sigma) - height);
				// temp [i * m + j] = (double) (rand () % 1000 - 500) / 100000.0 * (-vertical_grid [j] * vertical_grid [j] + config ["grid.z.width") * config ["grid.z.width") / 4.0);
			}
		}

		std::string file_format = "../input/" + config.get <std::string> ("input.file");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);

		io::output output_stream (new io::two_d::netcdf (n, m), buffer);

		double duration = 0.0;
		int mode = mode_flag;
		output_stream.append <double> ("T", &temp [0]);
		output_stream.append_scalar <double> ("t", &duration);
		output_stream.append_scalar <int> ("mode", &mode);

		output_stream.to_file ();
	} catch (std::exception& except) {
		std::cout << except.what () << '\n';
	}

	
	return 0;
}
