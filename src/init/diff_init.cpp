/*!***********************************************************************
 * \file init.cpp
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

/*!*******************************************************************
 * \brief The main call
 * 
 * \param argc The integer number of command line arguments
 * \param argv The character array of command line arguments
 *********************************************************************/
int main (int argc, char *argv[])
{
	try {
		int id = 0, n_elements = 1;

		// Initialize messenger
		mpi::messenger process_messenger (&argc, &argv);

		id = process_messenger.get_id ();
		n_elements = process_messenger.get_np ();
		
		io::parameters parameters = config (&argc, &argv, id);

		int m = parameters.get<int> ("grid.z.points") / n_elements + 1;
		m += (m - 1) % 2;
		int n = parameters.get<int> ("grid.x.points");
		double position_m0 = -parameters.get<double>("grid.z.width") / 2.0 + parameters.get<double>("grid.z.width") / n_elements * id;
		double position_mm = -parameters.get<double>("grid.z.width") / 2.0 + parameters.get<double>("grid.z.width") / n_elements * (id + 1);
		double position_n0 = -parameters.get<double>("grid.x.width") / 2.0;
		double position_nn = parameters.get<double>("grid.x.width") / 2.0;

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

		horizontal::grid horizontal_grid(new grids::axis(n, position_n0, position_nn));
		vertical::grid vertical_grid(new grids::axis(m, position_m0, position_mm, excess_0, excess_n));

		std::vector<double> temp_vec(n * m, 0.0), tempv_vec(n * m, 0.0);
		double* temp = &temp_vec[0];
		double* tempv = &tempv_vec[0];
		const double* pos_z = &vertical_grid[0];
		const double* pos_x = &horizontal_grid[0];

		double height = parameters.get<double>("grid.z.width") / 2.0;
		double sigma = 0.1;
		double scale = 1.0;
		double arg;
		// #pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				arg = (pos_z[j] * pos_z[j] + pos_x[i] * pos_x[i]) / 2. / sigma / sigma;
				temp[i * m + j] = scale * std::exp(-arg);
				tempv[i * m + j] = scale * (-pos_z[j] * pos_z[j] + height * height);
			}
		}

		std::string file_format = parameters.get<std::string> ("root") + parameters.get<std::string> ("input.directory") + parameters.get<std::string> ("input.file");
		char buffer[file_format.size() * 2];
		snprintf(buffer, file_format.size() * 2, file_format.c_str(), id);

		io::formatted_output<formats::netcdf> output_stream(formats::data_grid::two_d(n, m), buffer, formats::replace_file);

		double duration = 0.0;
		int mode = mode_flag;
		output_stream.append<double>("scalar", &temp[0]);
		output_stream.append<double>("x_velocity", &tempv[0]);

		output_stream.append<double>("t", &duration, formats::scalar);
		output_stream.append<int>("mode", &mode, formats::scalar);

		output_stream.to_file();
	} catch (std::exception& except) {
		FATAL (except.what());
	}

	
	return 0;
}
