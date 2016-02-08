/*!***********************************************************************
 * \file var_init.cpp
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

		int m = parameters.get <int> ("grid.z.points") / n_elements + 1;
		m += (m - 1) % 2;
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

		horizontal::grid horizontal_grid (new grids::axis (n, position_n0, position_nn));
		vertical::grid vertical_grid (new grids::axis (m, position_m0, position_mm, excess_0, excess_n));

		DEBUG ("TOTAL M: " << m);

		std::vector <double> temps_vec (n * m, 0.0), tempt_vec (n * m, 0.0);
		double *temps = &temps_vec [0], *tempt = &tempt_vec [0];
		const double *pos_z = &vertical_grid [0];
		
		double stop = 0.0, sbot = 1.0;
		stop = parameters.get <double> ("equations.composition.top.value", 0.0);
		sbot = parameters.get <double> ("equations.composition.bottom.value", 0.0);
		DEBUG ("stop " << stop << " sbot " << sbot);

		double ttop = 0.0, tbot = 1.0;
		// ttop = parameters.get <double> ("equations.temperature.top.value", 0.0);
		tbot = parameters.get <double> ("equations.temperature.bottom.value", 0.0);
		ttop = tbot + parameters.get <double> ("equations.temperature.top.value", 0.0) * parameters.get <double> ("grid.z.width");

		double height = parameters.get <double> ("grid.z.width");
		double diff_bottom = parameters.get <double> ("equations.temperature.diffusion");
		double diff = parameters.get <double> ("equations.temperature.diffusion");
		double mu_length = parameters.get <double> ("equations.temperature.tanh.length");
		double mu_zero = parameters.get <double> ("equations.temperature.tanh.zero");
		double stiffness = parameters.get <double> ("equations.temperature.stiffness");
		double chi = 1. / (1. - parameters.get <double> ("equations.temperature.sources.z_velocity"));

		if (chi < 0.) {
			FATAL("Chi cannot be less than zero.");
			throw 3;
		}

		double a = chi - 1.;
		double b = 1. + 2. * (stiffness - 2.) * chi + chi * chi;
		double c = 2. * (1. + stiffness) * chi;
		double d = chi * (1. - stiffness);
		double e = chi * (1 + stiffness);
		
		double scale = 0.001;
		double width = parameters.get <double> ("grid.z.width");
		double arg;


		INFO("Chi = " << chi);
		INFO("Stiffness = " << stiffness);
		INFO("Diffusion = " << diff);
		INFO("");

		INFO("a = " << a);
		INFO("b = " << b);
		INFO("c = " << c);
		INFO("d = " << d);
		INFO("e = " << e);
		// if (b + c < 0.) {
		// 	FATAL("Problem setup failed at lower boundary: try different values of chi, stiffness, and diffusion");
		// 	throw 1;
		// }
		// if (b - c < 0.) {
		// 	FATAL("Problem setup failed at upper boundary: try different values of chi, stiffness, and diffusion");
		// 	throw 2;
		// }

		double phi = ((1.0 - chi - sqrt(1.0 - 2.0 * chi - 4.0 * stiffness * chi + chi * chi)) / (2.0 * stiffness * chi) - 1.0) / 2.0;
		if (phi < 0.0) {
			FATAL("Cannot realize stiffness with positive diffusion. Try changing chi or the stiffness.")
			throw 2;
		}
		WARN("Phi = " << phi);

		// #pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				// double temp = (-std::cosh((pos_z[j] + mu_zero) / mu_length) + phi*std::sinh((pos_z[j] + mu_zero) / mu_length)) / (-std::cosh((-1. + mu_zero) / mu_length) + phi*std::sinh((-1. + mu_zero) / mu_length));
				// tempt [i * m + j] = (-tbot*chi + tbot*phi*phi*chi + (1. - phi*std::tanh((1. + mu_zero)/mu_length))*(1. + pos_z[j] + mu_length*phi*log(temp))) / (-1. + phi*phi) / chi;
				
				// tempt [i * m + j] = tbot + (pos_z[j] + height / 2) * (ttop - tbot) / height;
				// 
				if (j == 0) {
					tempt[i * m + j] = tbot;
				} else {
					arg = (-(pos_z[j] + pos_z[j - 1]) / 2.0 - mu_zero) / mu_length;
					tempt[i * m + j] = tempt[i * m + j - 1] - (pos_z[j] - pos_z[j - 1]) * (1. - phi) / chi / (1.0 + phi * std::tanh(arg));
					// INFO ((-diff / (diff + (a - sqrt(b - c * std::tanh(arg))) / (d + e * std::tanh(arg))) * (chi - 1.)) / chi << " " << (b - c * std::tanh(arg)));
				}
			}
		}

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				temps [i * m + j] = (stop - sbot) / (height) * (pos_z [j] + height / 2.0) + sbot;

				temps [i * m + j] += (double) (rand () % 2000 - 1000) * scale / 1.0e3 * std::cos(pos_z [j] * 3.14159 / width);
				// tempt [i * m + j] += (double) (rand () % 2000 - 1000) * scale / 1.0e3 * std::cos(pos_z [j] * 3.14159 / width);
			}
		}

		std::string file_format = parameters.get <std::string> ("root") + parameters.get <std::string> ("input.directory") + parameters.get <std::string> ("input.file");
		char buffer [file_format.size () * 2];
		snprintf (buffer, file_format.size () * 2, file_format.c_str (), id);

		io::formatted_output <formats::netcdf> output_stream (formats::data_grid::two_d (n, m), buffer, formats::replace_file);

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
