/*!***********************************************************************
 * \file element_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_3SURDTOH
#define ELEMENT_HPP_3SURDTOH

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "../config.hpp"
#include "../bases/element.hpp"
#include "../bases/plan.hpp"
#include "../utils/utils.hpp"
#include "../utils/chebyshev.hpp"
#include "../utils/interpolate.hpp"
	
namespace one_d
{
	enum edges {
		edge_0 = 0,
		edge_n = 1
	};
	
	/*!*******************************************************************
	 * A 1D implementation of the element base class. This provides the
	 * storage, indexing facilities, and failsafe_dump output. The plans should be added in a 
	 * further subclass.
	 * 
	 * \brief \copybrief bases::element
	 *********************************************************************/
	class element : public bases::element
	{
	public:	
		/*!*******************************************************************
		 * \param i_n The number of data elements in each scalar
		 * \param i_position_0 The double position of index excess_0
		 * \param i_position_n The double position of index n - 1 - excess_n
		 * \param i_excess_0 The integer number of points evaluated in the adjacent element
		 * \param i_excess_n The integer number of points evaluated in the adjacent element
		 * \copydoc bases::element::element ()
		 *********************************************************************/
		element (int i_n, double i_position_0, double i_position_n, int i_name, io::parameter_map& i_inputParams, utils::messenger* i_messenger_ptr, int i_flags) : 
		bases::element (i_name, 2, i_inputParams, i_messenger_ptr, i_flags) {
			n = i_n;
			position_0 = i_position_0;
			position_n = i_position_n;
			boundary_weights [edge_0] = 0.0;
			boundary_weights [edge_n] = 0.0;
			
			cell.resize (i_n);
			for (int i = 0; i < i_n; ++i) {
				cell [i] = i;
			}
			
			failsafe_dump = std::make_shared <io::simple_output> (io::simple_output ("dump_" + std::to_string (name) + ".dat", n));
			failsafe_dump->append (&cell [0]);
		}
		
		virtual ~element () {}
	
		/*!*******************************************************************
		 * \copydoc bases::element::initialize ()
		 *********************************************************************/
		virtual void initialize (int name, double* initial_conditions = NULL) {
			if (scalars [name].size () == (unsigned int) 0) {
				names.push_back (name);
				scalars [name].resize (n, 0.0);
			}
			if (initial_conditions) {
				utils::copy (n, initial_conditions, &(scalars [name]) [0]);
			}
			fixed_points_0 [name] = scalars [name] [0];
			fixed_points_n [name] = scalars [name] [n - 1];
			failsafe_dump->append (&(scalars [name]) [0]);
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element::operator[] ()
		 *********************************************************************/
		inline double& operator[] (int name) {
			if (scalars [name].size () == (unsigned int) 0) {
				initialize (name);
			}
			return scalars [name] [0];
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element::explicit_reset ();
			bases::element::iterator iter;
			for (iter = begin (); iter != end (); ++iter) {
				if (*iter < 0) {
					utils::scale (n, 0.0, &((*this) [*iter]));
				}
			}
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element::execute_boundaries ()
		 *********************************************************************/
		inline void execute_boundaries () {
			if (!(boundary_bools [edge_0])) {
				for (iterator iter = begin (); iter != end (); ++iter) {
					(*this) (*iter, 0) = fixed_points_0 [*iter];
				}
			}
			if (!(boundary_bools [edge_n])) {
				for (iterator iter = begin (); iter != end (); ++iter) {
					(*this) (*iter, n - 1) = fixed_points_n [*iter];
				}
			}
		}
		
	protected:
		int n; //!< The number of elements in each 1D array
		double position_0; //!< The double position of index 0
		double position_n; //!< The double position of index n - 1
		std::vector<int> cell; //!< An integer array for tracking each cell number for output

		std::map <int, std::vector <double>> scalars; //!< A vector of scalar vectors
		std::map <int, double> fixed_points_0; //!< The initial values of the scalars at index 0
		std::map <int, double> fixed_points_n; //!< The initial values of the scalars at index n - 1
		
		/*
			TODO Perhaps there's a better way to handle the fixed points
		*/
	};

	namespace chebyshev
	{
		/*!*******************************************************************
		 * \brief A Chebyshev implementation of the 1D element class
		 *********************************************************************/
		class element : public one_d::element
		{
		public:
			/*!*******************************************************************
			 * \copydoc one_d::element::element ()
			 *********************************************************************/
			element (int i_n, double i_position_0, double i_position_n, int i_name, io::parameter_map& i_inputParams, utils::messenger* i_messenger_ptr, int i_flags) : one_d::element (i_n, i_position_0, i_position_n, i_name, i_inputParams, i_messenger_ptr, i_flags) {
				initialize (position);
				set_grid (std::make_shared<chebyshev_grid> (chebyshev_grid (i_n, i_n, sqrt (2.0 / (i_n - 1.0)), position_0 - position_n)));
			}
			virtual ~element () {}
				
			/*!*******************************************************************
			 * \copydoc one_d::element::initialize ()
			 *********************************************************************/
			virtual void initialize (int name, double* initial_conditions = NULL) {
				TRACE ("Initializing " << name);
				if (name == position && !initial_conditions) {
					double pioN = std::acos (-1.0) / (n - 1);
					double scale = (position_0 - position_n) / 2.0;
					double initial_position = (position_0 + position_n) / 2.0;
					std::vector <double> init (n);
					for (int i = 0; i < n; ++i) {
						init [i] = scale * std::cos (i * pioN) + initial_position;
					}
					one_d::element::initialize (name, &init [0]);
				} else if (name == velocity && !initial_conditions){
					double scale = inputParams["init_cond_scale"].asDouble;
					double width = inputParams["init_cond_width"].asDouble;
					double mean = inputParams["init_cond_mean"].asDouble;
					double sigma = inputParams["init_cond_sigma"].asDouble;
					std::vector <double> init (n);
					double height, temp;
					height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
					for (int i = 0; i < n; ++i) {
						temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
						if (temp > 0.0) {
							init [i] = temp;
						} else {
							init [i] = 0.0;
						}
					}
					one_d::element::initialize (name, &init [0]);
				} else {
					one_d::element::initialize (name, initial_conditions);
				}
				TRACE ("Initialized.");
			}
		};
		
		/*!*******************************************************************
		 * \brief A simple implementation of the element class with diffusion
		 * 
		 * This class contains a full element's capacity to run a single 
		 * element diffusion in 1D with constant timestep.
		 *********************************************************************/
		class advection_diffusion_element : public element
		{
		public:
			/*!*******************************************************************
			 * \copydoc element::element ()
			 *********************************************************************/
			advection_diffusion_element (int i_n, double i_position_0, double i_position_n, int i_excess_0, int i_excess_n, int i_name, io::parameter_map& i_inputParams, utils::messenger* i_messenger_ptr, int i_flags);
			
			virtual ~advection_diffusion_element () {}
		
			/*!*******************************************************************
			 * \copydoc element::implicit_reset ()
			 *********************************************************************/
			inline void implicit_reset () {
				element::implicit_reset ();
			
				if (!(flags & factorized)) {
					utils::scale (n * n, 0.0, &matrix [0]);
				}
			}
			
			virtual double calculate_timestep ();
		
		private:
			std::vector<double> matrix; //!< A vector containing the double matrix used in the implicit solver
		};
		
		class cuda_element : public element
		{
		public:
			cuda_element (int i_n, double i_position_0, double i_position_n, int i_excess_0, int i_excess_n, int i_name, io::parameter_map& i_input_Params, utils::messenger* i_messenger_ptr, int i_flags);
			
			virtual ~cuda_element () {}
		
			inline void implicit_reset () {
				if (!(flags & factorized)) {
					utils::scale (n * n, 0.0, &matrix [0]);
				}
			}
		
		private:
			std::vector<double> matrix; //!< A vector containing the double matrix used in the implicit solver
		};
	} /* chebyshev */
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
