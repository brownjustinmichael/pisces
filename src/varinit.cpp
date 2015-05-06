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
#include "plans/grids/grid.hpp"
#include "config.hpp"

using namespace grids;

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

		horizontal::grid <double> horizontal_grid (new grids::axis (n, position_n0, position_nn));
		vertical::grid <double> vertical_grid (new grids::axis (m, position_m0, position_mm, excess_0, excess_n));

		std::vector <double> temps_vec (n * m, 0.0), tempt_vec (n * m, 0.0);
		double *temps = &temps_vec [0], *tempt = &tempt_vec [0];
		const double *pos_z = &vertical_grid [0];
		
		double stop = 0.0, sbot = 0.0;
		if (parameters ["equations.composition.top.type"].as <std::string> () == "fixed_value") {
			stop = parameters ["equations.composition.top.value"].as <double> ();
		}
		if (parameters ["equations.composition.bottom.type"].as <std::string> () == "fixed_value") {
			sbot = parameters ["equations.composition.bottom.value"].as <double> ();
		}
		
		double ttop = 0.0, tbot = 0.0;
		if (parameters ["equations.temperature.top.type"].as <std::string> () == "fixed_value") {
			ttop = parameters ["equations.temperature.top.value"].as <double> ();
		}
		if (parameters ["equations.temperature.bottom.type"].as <std::string> () == "fixed_value") {
			tbot = parameters ["equations.temperature.bottom.value"].as <double> ();
		}
		
		double height = parameters.get <double> ("grid.z.width");
		double diff_bottom = parameters.get <double> ("equations.temperature.diffusion");
		double diff_top = parameters ["equations.temperature.bg_diffusion"].IsDefined () ? parameters.get <double> ("equations.temperature.bg_diffusion") : diff_bottom;
		double diff_width = parameters ["equations.temperature.diff_width"].IsDefined () ? parameters.get <double> ("equations.temperature.diff_width") : parameters.get <double> ("equations.composition.diff_width");
		double bg_scale = -(tbot - ttop) / ((height / 2. - diff_width) * (1. / diff_top + 1. / diff_bottom) + 2. * diff_width / (diff_top - diff_bottom) * log (diff_top / diff_bottom));
		
		double scale = parameters.get <double> ("init.scale");
		#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				temps [i * m + j] = sbot;
				// tempt [i * m + j] = (ttop - tbot) / (height) * (pos_z [j] + height / 2.0) + tbot;
				tempt [i * m + j] = tbot + bg_scale * (std::min (pos_z [j], -diff_width) + height / 2.) / diff_bottom;
				if (pos_z [j] > -diff_width) {
					tempt [i * m + j] += bg_scale * 2. * diff_width / (diff_top - diff_bottom) * log ((diff_top - diff_bottom) / 2. / diff_width / diff_bottom * std::min (pos_z [j], diff_width) + (diff_top + diff_bottom) / 2. / diff_bottom);
					temps [i * m + j] = (stop - sbot) / (diff_width) * (pos_z [j]) / 2.0;
					if (pos_z [j] > diff_width) {
						tempt [i * m + j] += bg_scale * (pos_z [j] - diff_width) / diff_top;
						temps [i * m + j] = stop;
					}
				}
				temps [i * m + j] += (double) (rand () % 2000 - 1000) * scale / 1.0e3;
				tempt [i * m + j] += (double) (rand () % 2000 - 1000) * scale / 1.0e3;
			}
		}

		std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + parameters.get <std::string> ("input.file");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);

		io::formatted_output <formats::netcdf> output_stream (formats::data_grid::two_d (n, m, 0, parameters.get <bool> ("input.full") ? n_elements * m : 0, 0, parameters.get <bool> ("input.full") ? id * m : 0), buffer, formats::replace_file);

		double duration = 0.0;
		int mode = mode_flag;
		output_stream.append <double> ("temperature", &tempt [0]);
		output_stream.append <double> ("composition", &temps [0]);

		output_stream.append <double> ("t", &duration, formats::scalar);
		output_stream.append <int> ("mode", &mode, formats::scalar);

		output_stream.to_file ();
	} catch (std::exception& except) {
		FATAL (except.what ());
	}

	
	return 0;
}
