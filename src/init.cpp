/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <sstream>
#include <iomanip>

#include "logger/logger.hpp"
#include "io/output.hpp"
#include "io/parameters.hpp"
#include "io/formats/netcdf.hpp"
#include "mpi/messenger.hpp"
#include "plans/grid.hpp"

using namespace plans;

int main (int argc, char *argv[])
{
	try {
		int id = 0, n_elements = 1;

		// Initialize messenger
		mpi::messenger process_messenger (&argc, &argv);

		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();
		
		logger::log_config::configure (&argc, &argv, id, "init_%d.log");

		std::string config_filename;

		if (argc <= 1) {
			config_filename = "config.yaml";
		} else {
			config_filename = argv [1];
		}
		io::parameters config (config_filename);
		if (!config ["grid.x.points"].IsDefined ()) config ["grid.x.points"] = 64;
		if (!config ["grid.z.points"].IsDefined ()) config ["grid.z.points"] = 64;

		int m = config.get <int> ("grid.z.points") / n_elements + 1;
		m += m % 2;
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

		horizontal::grid <double> horizontal_grid (new plans::axis (n, position_n0, position_nn));
		vertical::grid <double> vertical_grid (new plans::axis (m, position_m0, position_mm, excess_0, excess_n));

		std::vector <double> temps (n * m, 0.0), tempt (n * m, 0.0);

		double scale = config.get <double> ("init.scale");
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				// temps [i * m + j] = scale * (vertical_grid [j]);
				temps [i * m + j] += (double) (rand () % 1000 - 500) / 100000.0 * (-vertical_grid [j] * vertical_grid [j] + config.get <double> ("grid.z.width") * config.get <double> ("grid.z.width") / 4.0) / config.get <double> ("grid.z.width") / config.get <double> ("grid.z.width");
				tempt [i * m + j] += (double) (rand () % 1000 - 500) / 100000.0 * (-vertical_grid [j] * vertical_grid [j] + config.get <double> ("grid.z.width") * config.get <double> ("grid.z.width") / 4.0) / config.get <double> ("grid.z.width") / config.get <double> ("grid.z.width");
				tempt [i * m + j] += scale * sin (8 * 3.14159 / config.get <double> ("grid.z.width") * vertical_grid [j]) * sin (4 * 3.14159 / config.get <double> ("grid.x.width") * horizontal_grid [i]);
				// temp [i * m + j] = 0.01 * std::cos (std::acos (-1.0) * vertical_grid [j]) * std::sin (8.0 * std::acos (-1.0) * horizontal_grid [i] / 3.0);
			}
		}

		std::string file_format = "input/" + config.get <std::string> ("input.file");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);

		io::formatted_output <io::formats::two_d::netcdf> output_stream (io::data_grid::two_d (n, m, 0, config.get <bool> ("input.full") ? n_elements * m : 0, 0, config.get <bool> ("input.full") ? id * m : 0), buffer, io::replace_file);

		double duration = 0.0;
		int mode = mode_flag;
		output_stream.append <double> ("T", &tempt [0]);
		output_stream.append <double> ("S", &temps [0]);

		output_stream.append <double> ("t", &duration, io::scalar);
		output_stream.append <int> ("mode", &mode, io::scalar);

		output_stream.to_file ();
	} catch (std::exception& except) {
		std::cout << except.what () << '\n';
	}

	
	return 0;
}
