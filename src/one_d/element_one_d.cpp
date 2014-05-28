/*!***********************************************************************
 * \file element_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "element_one_d.hpp"

#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <iomanip>
#include "../config.hpp"
#include "diffusion_one_d.hpp"
#include "advection_one_d.hpp"
#include "solver_one_d.hpp"
#include "transform_one_d.hpp"
	
namespace one_d
{
	namespace cosine
	{
		template <class datatype>
		int element <datatype>::mode = mode_flag;
		
		template <class datatype>
		advection_diffusion_element <datatype>::advection_diffusion_element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
		element <datatype> (i_axis_n, i_name, i_params, i_messenger_ptr, i_element_flags) {
		
			assert (n > 0);
		
			TRACE ("Initializing...");
			io::parameters& nparams = params;
			advection_coeff = nparams.get <datatype> ("velocity.advection");
			cfl = nparams.get <datatype> ("time.cfl");
			datatype scale = nparams.get <datatype> ("init.scale");
			datatype mean = nparams.get <datatype> ("init.mean");
			datatype sigma = nparams.get <datatype> ("init.sigma");
			datatype width = nparams.get <datatype> ("grid.z.width");
			std::vector <datatype> init (n);
			datatype height, temp;
			height = std::max (scale * std::exp (- (-width / 2.0 - mean) * (-width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] = temp;
				} else {
					init [i] = 0.0;
				}
			}
			position_ptr = ptr (position);
			velocity_ptr = initialize (velocity, "u", &init [0]);

			// // Set up output
			// std::string file_format = "../output/" + i_params.get <std::string> ("output.file");
			// char buffer [file_format.size () * 2];
			// snprintf (buffer, file_format.size () * 2, file_format.c_str (), name);
			// 
			// normal_stream.reset (new io::incremental (new io::one_d::netcdf (n), buffer, i_params.get <int> ("output.every")));
			// normal_stream->template append <int> ("i", &(cell [0]));
			// normal_stream->template append <datatype> ("x", ptr (position));
			// normal_stream->template append <datatype> ("u", ptr (velocity));
			
			// // Set up solver
			element <datatype>::add_solver (velocity, new solver <datatype> (*grids [0], messenger_ptr, timestep, alpha_0, alpha_n, ptr (velocity), &element_flags [state], &element_flags [velocity]));
			// 
			// // Set up plans in order
			solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new diffusion <datatype> (*solvers [velocity], nparams.get <datatype> ("velocity.diffusion"), nparams.get <datatype> ("time.alpha"))), pre_plan);
			if (nparams.get <datatype> ("velocity.advection") != 0.0) {
				solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new advection <datatype> (*solvers [velocity], nparams.get <datatype> ("velocity.advection"))), post_plan);
			}
								
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype advection_diffusion_element <datatype>::calculate_timestep (int i, io::virtual_dump *dump) {
			return (std::abs ((position_ptr [i - 1] - position_ptr [i + 1]) / velocity_ptr [i] / advection_coeff) * cfl);
		}
		
		template class advection_diffusion_element <double>;
	}
	
	namespace chebyshev
	{
		template <class datatype>
		int element <datatype>::mode = mode_flag;
		
		template <class datatype>
		advection_diffusion_element <datatype>::advection_diffusion_element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
		element <datatype> (i_axis_n, i_name, i_params, i_messenger_ptr, i_element_flags) {
		
			assert (n > 0);
			io::parameters& nparams = params;
			advection_coeff = nparams.get <datatype> ("velocity.advection");
			cfl = nparams.get <datatype> ("time.cfl");
		
			TRACE ("Initializing...");
			datatype scale = nparams.get <datatype> ("init.scale");
			datatype mean = nparams.get <datatype> ("init.mean");
			datatype sigma = nparams.get <datatype> ("init.sigma");
			datatype width = nparams.get <datatype> ("grid.z.width");
			std::vector <datatype> init (n);
			datatype height, temp;
			height = std::max (scale * std::exp (- (-width / 2.0 - mean) * (-width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] = temp;
				} else {
					init [i] = 0.0;
				}
			}
			position_ptr = ptr (position);
			velocity_ptr = initialize (velocity, "u", &init [0]);
			
			// // Set up output
			// std::ostringstream filestream;
			// filestream << "../output/" + nparams.get <std::string> ("output.file") << "_" << std::setfill ('0') << std::setw (2) << name << "_%04i";
			// normal_stream.reset (new io::incremental (new io::one_d::netcdf (n), filestream.str (), nparams.get <int> ("output.every")));
			// normal_stream->template append <int> ("i", &(cell [0]));
			// normal_stream->template append <datatype> ("x", ptr (position));
			// normal_stream->template append <datatype> ("u", ptr (velocity));
			
			// // Set up solver
			element <datatype>::add_solver (velocity, new solver <datatype> (*grids [0], messenger_ptr, timestep, alpha_0, alpha_n, ptr (velocity), &element_flags [state], &element_flags [velocity]));
			// 
			// // Set up plans in order
			solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new diffusion <datatype> (*solvers [velocity], nparams.get <datatype> ("velocity.diffusion"), nparams.get <datatype> ("time.alpha"))), pre_plan);
			if (nparams.get <datatype> ("velocity.advection") != 0.0) {
				solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new advection <datatype> (*solvers [velocity], nparams.get <datatype> ("velocity.advection"))), post_plan);
			}
			
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype advection_diffusion_element <datatype>::calculate_timestep (int i, io::virtual_dump *dump) {
			return (std::abs ((position_ptr [i - 1] - position_ptr [i + 1]) / velocity_ptr [i] / advection_coeff) * cfl);
		}
		
		template class advection_diffusion_element <double>;
		
		template <class datatype>
		nonlinear_diffusion_element <datatype>::nonlinear_diffusion_element (bases::axis i_axis_n, int i_name, io::parameters& params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
		element <datatype> (i_axis_n, i_name, params, i_messenger_ptr, i_element_flags) {
			io::parameters& nparams = params;
			advection_coeff = nparams.get <datatype> ("velocity.advection");
			cfl = nparams.get <datatype> ("time.cfl");
			datatype diffusion_coeff = nparams.get <datatype> ("velocity.diffusion");
			datatype alpha = nparams.get <datatype> ("time.alpha");
		
			assert (n > 0);
		
			TRACE ("Initializing...");
			
			datatype scale = nparams.get <datatype> ("init.scale");
			datatype width = nparams.get <datatype> ("grid.width");
			datatype mean = nparams.get <datatype> ("init.mean") - 0.5;
			datatype sigma = nparams.get <datatype> ("init.sigma");
			std::vector <datatype> init (n);
			datatype height, temp;
			for (int i = 0; i < n; ++i) {
				// if ((*this) (position, i) < 0.0) {
				// 	init [i] = 0.0;
				// } else {
				// 	init [i] = 1.0;
				// }
				init [i] = ((*this) (position, i) - (*this) (position)) * scale;
			}
			height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] += temp;
				}
			}
			mean = nparams.get <datatype> ("init.mean") + 0.5;
			height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
			for (int i = 0; i < n; ++i) {
				temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
				if (temp > 0.0) {
					init [i] += temp;
				}
			}
			position_ptr = ptr (position);
			velocity_ptr = initialize (velocity, "u", &init [0]);
			
			// // Set up output
			// std::ostringstream filestream;
			// filestream << "../output/" + nparams.get <std::string> ("output.file") + "_" << std::setfill ('0') << std::setw (2) << name << "_%04i";
			// normal_stream.reset (new io::incremental (new io::one_d::netcdf (n), filestream.str (), nparams.get <int> ("output.every")));
			// normal_stream->template append <int> ("i", &(cell [0]));
			// normal_stream->template append <datatype> ("x", ptr (position));
			// normal_stream->template append <datatype> ("u", ptr (velocity));			

			// Set up solver
			element <datatype>::add_solver (velocity, new solver <datatype> (*grids [0], messenger_ptr, timestep, alpha_0, alpha_n, ptr (velocity), &element_flags [state], &element_flags [velocity]));
			
			// Set up plans in order
			solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new diffusion <datatype> (*solvers [velocity], diffusion_coeff, alpha)), pre_plan);
			if (nparams.get <datatype> ("velocity.nonlinear") != 0.0) {
				solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new nonlinear_diffusion <datatype> (*solvers [velocity], nparams.get <datatype> ("velocity.nonlinear"))), post_plan);
			}
			if (advection_coeff != 0.0) {
				solvers [velocity]->add_plan (std::shared_ptr <bases::plan <datatype>> (new advection <datatype> (*solvers [velocity], advection_coeff)), post_plan);
			}
		
			// normal_stream->to_file ();
		
			TRACE ("Initialized.");
		}
		
		template <class datatype>
		datatype nonlinear_diffusion_element <datatype>::calculate_timestep (int i, io::virtual_dump *dump) {
			io::parameters& nparams = params;
			datatype advection_coeff = nparams.get <datatype> ("velocity.advection");
			datatype cfl = nparams.get <datatype> ("time.cfl");
			datatype nonlinear = nparams.get <datatype> ("velocity.nonlinear");
			datatype t_timestep = std::min (nparams.get <datatype> ("time.max"), std::abs ((position_ptr [i - 1] - position_ptr [i + 1]) / velocity_ptr [i] / advection_coeff) * cfl);
			return std::min (t_timestep, std::abs ((*this) (position, i + 1) - (*this) (position, i - 1)) * ((*this) (position, i + 1) - (*this) (position, i - 1)) / 2.0 / nonlinear / (*this) (velocity, i) * cfl);
		}
		
		template class nonlinear_diffusion_element <double>;
	} /* chebyshev */
} /* one_d */