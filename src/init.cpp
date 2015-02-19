/*!***********************************************************************
 * \file main.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "mpi/messenger.hpp"

#include <sstream>
#include <iomanip>

#include "logger/logger.hpp"
#include "io/output.hpp"
#include "io/parameters.hpp"
#include "io/formats/netcdf.hpp"
#include "plans/grid.hpp"
#include "config.hpp"

using namespace plans;

int main (int argc, char *argv[])
{
	try {
		int id = 0, n_elements = 1;

		// Initialize messenger
		mpi::messenger process_messenger (&argc, &argv);

		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();
		
		io::parameters parameters = config (&argc, &argv, id);

		int m = parameters.get <int> ("grid.z.points") / n_elements + 1;
		m += m % 2;
		int n = parameters.get <int> ("grid.x.points");
		double position_m0 = -parameters.get <double> ("grid.z.width") / 2.0 + parameters.get <double> ("grid.z.width") / n_elements * id;
		double position_mm = -parameters.get <double> ("grid.z.width") / 2.0 + parameters.get <double> ("grid.z.width") / n_elements * (id + 1);
		double position_n0 = -parameters.get <double> ("grid.x.width") / 2.0;
		double position_nn = parameters.get <double> ("grid.x.width") / 2.0;

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

		std::vector <double> temps_vec (n * m, 0.0), tempt_vec (n * m, 0.0);
		double *temps = &temps_vec [0], *tempt = &tempt_vec [0];

		double scale = parameters.get <double> ("init.scale");
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				// temps [i * m + j] = scale * (vertical_grid [j]);
				temps [i * m + j] += (double) (rand () % 1000 - 500) / 100000.0 * (-vertical_grid [j] * vertical_grid [j] + parameters.get <double> ("grid.z.width") * parameters.get <double> ("grid.z.width") / 4.0) / parameters.get <double> ("grid.z.width") / parameters.get <double> ("grid.z.width");
				tempt [i * m + j] += (double) (rand () % 1000 - 500) / 100000.0 * (-vertical_grid [j] * vertical_grid [j] + parameters.get <double> ("grid.z.width") * parameters.get <double> ("grid.z.width") / 4.0) / parameters.get <double> ("grid.z.width") / parameters.get <double> ("grid.z.width");
				tempt [i * m + j] += scale * sin (8 * 3.14159 / parameters.get <double> ("grid.z.width") * vertical_grid [j]) * sin (4 * 3.14159 / parameters.get <double> ("grid.x.width") * horizontal_grid [i]);
				temps [i * m + j] += scale * sin (8 * 3.14159 / parameters.get <double> ("grid.z.width") * vertical_grid [j]) * sin (4 * 3.14159 / parameters.get <double> ("grid.x.width") * horizontal_grid [i]);
				// temp [i * m + j] = 0.01 * std::cos (std::acos (-1.0) * vertical_grid [j]) * std::sin (8.0 * std::acos (-1.0) * horizontal_grid [i] / 3.0);
			}
		}

		std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + parameters.get <std::string> ("input.file");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);

		io::formatted_output <io::formats::two_d::netcdf> output_stream (io::data_grid::two_d (n, m, 0, parameters.get <bool> ("input.full") ? n_elements * m : 0, 0, parameters.get <bool> ("input.full") ? id * m : 0), buffer, io::replace_file);

		double duration = 0.0;
		int mode = mode_flag;
		output_stream.append <double> ("temperature", &tempt [0]);
		output_stream.append <double> ("composition", &temps [0]);

		output_stream.append <double> ("t", &duration, io::scalar);
		output_stream.append <int> ("mode", &mode, io::scalar);

		output_stream.to_file ();
	} catch (std::exception& except) {
		FATAL (except.what ());
	}

	
	return 0;
}
