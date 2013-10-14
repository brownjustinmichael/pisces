/*!***********************************************************************
 * \file element_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_3SURDTOH
#define ELEMENT_HPP_3SURDTOH

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include "../config.hpp"
#include "../bases/element.hpp"
#include "../utils/utils.hpp"
#include "../bases/grid.hpp"
	
namespace one_d
{
	/*!**********************************************************************
	 * \brief Integer representation of the edges
	 * 
	 * These edges must range from 0 to the total number of boundaries - 1 for
	 * messenger to work correctly.
	 ************************************************************************/
	enum edges {
		edge_0 = 0,
		edge_n = 1
	};
	
	/*!*******************************************************************
	 * \brief The 1D base element class
	 * 
	 * A 1D implementation of the element base class. This provides the
	 * storage, indexing facilities, and failsafe_dump output. The plans should be added in a 
	 * further subclass.
	 *********************************************************************/
	template <class datatype>
	class element : public bases::element <datatype>
	{
	public:	
		/*!*******************************************************************
		 * \param i_n The number of data elements in each scalar
		 * \param i_position_0 The datatype position of index excess_0
		 * \param i_position_n The datatype position of index n - 1 - excess_n
		 * \copydoc bases::element <datatype>::element ()
		 *********************************************************************/
		element (bases::axis *i_axis_n, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 1, i_inputParams, i_messenger_ptr, i_flags),
		axis_n (i_axis_n),
		n (axis_n->n) {
			cell.resize (n);
			for (int i = 0; i < n; ++i) {
				cell [i] = i;
			}
			
			edge_map [edge_0] = 0;
			edge_map [edge_n] = n - 1;
			
			if (i_messenger_ptr->linked (edge_0)) {
				fixed_points [edge_0];
			}
			if (i_messenger_ptr->linked (edge_n)) {
				fixed_points [edge_n];
			}
			
			std::ostringstream convert;
			convert << name;
			failsafe_dump.reset (new io::simple_output <datatype>  ("dump_" + convert.str () + ".dat", n));
			failsafe_dump->append (&cell [0]);
		}
		
		virtual ~element () {}
	
		/*!*******************************************************************
		 * \brief Get the datatype reference to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * \return A datatype reference to the first element of the named scalar
		 *********************************************************************/
		inline datatype& operator[] (int name) {
			if (scalars [name].size () == (unsigned int) 0) {
				initialize (name);
			}
			return scalars [name] [0];
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual void initialize (int name, datatype* initial_conditions = NULL) {
			if (scalars [name].size () == (unsigned int) 0) {
				scalars [name].resize (n, 0.0);
			}
			if (initial_conditions) {
				utils::copy (n, initial_conditions, &(scalars [name]) [0]);
			}
			for (std::map <int, int>::iterator i_edge = edge_map.begin (); i_edge != edge_map.end (); ++i_edge) {
				fixed_points [i_edge->first] [name] = scalars [name] [i_edge->second];
			}
			failsafe_dump->append (&(scalars [name]) [0]);
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element <datatype>::explicit_reset (); 
			for (typename std::map <int, std::vector <datatype> >::iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				if (iter->first < 0) {
					utils::scale (n, 0.0, &(iter->second [0]));
				}
			}
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::execute_boundaries ()
		 *********************************************************************/
		inline void execute_boundaries () {
			for (typename std::map <int, std::map <int, datatype> >::iterator i_edge = fixed_points.begin (); i_edge != fixed_points.end (); ++i_edge) {
				for (typename std::map <int, datatype>::iterator iter = (i_edge->second).begin (); iter != (i_edge->second).end (); ++iter) {
					(*this) (iter->first, edge_map [i_edge->first]) = iter->second;
				}
			}
		}
		
	protected:
		using bases::element <datatype>::name;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		
		bases::axis *axis_n;
		int &n;
		std::vector<int> cell; //!< An integer array for tracking each cell number for output

		std::map <int, std::vector <datatype>> scalars; //!< A vector of scalar vectors
		
		std::map <int, int> edge_map;
		std::map <int, std::map <int, datatype> > fixed_points;
	};

	namespace chebyshev
	{
		/*!*******************************************************************
		 * \brief A Chebyshev implementation of the 1D element class
		 *********************************************************************/
		template <class datatype>
		class element : public one_d::element <datatype>
		{
		public:
			/*!*******************************************************************
			 * \copydoc one_d::element::element ()
			 *********************************************************************/
			element (bases::axis *i_axis_n, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_axis_n, i_name, i_inputParams, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				one_d::element <datatype>::set_grid (new bases::chebyshev::grid <datatype> (axis_n, sqrt (2.0 / (n - 1.0)), messenger_ptr->linked (edge_0), messenger_ptr->linked (edge_n)));
				initialize (position);
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
				
			/*!*******************************************************************
			 * \copydoc one_d::element::initialize ()
			 *********************************************************************/
			virtual void initialize (int name, datatype* initial_conditions = NULL) {
				TRACE ("Initializing " << name);
				if (name == position && !initial_conditions) {
					one_d::element <datatype>::initialize (name, &(grids [0]->position ()));
				} else if (name == velocity && !initial_conditions){
					datatype scale = inputParams["init_cond_scale"].asDouble;
					datatype width = inputParams["init_cond_width"].asDouble;
					datatype mean = inputParams["init_cond_mean"].asDouble;
					datatype sigma = inputParams["init_cond_sigma"].asDouble;
					std::vector <datatype> init (n);
					datatype height, temp;
					height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
					for (int i = 0; i < n; ++i) {
						temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
						if (temp > 0.0) {
							init [i] = temp;
						} else {
							init [i] = 0.0;
						}
					}
					one_d::element <datatype>::initialize (name, &init [0]);
				} else {
					one_d::element <datatype>::initialize (name, initial_conditions);
				}
				TRACE ("Initialized.");
			}
			
		protected:
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::axis_n;
			using one_d::element <datatype>::grids;
			using one_d::element <datatype>::inputParams;
			using one_d::element <datatype>::messenger_ptr;
		};
		
		/*!*******************************************************************
		 * \brief A simple implementation of the element class with diffusion
		 * 
		 * This class contains a full element's capacity to run a single 
		 * element diffusion in 1D with constant timestep.
		 *********************************************************************/
		template <class datatype>
		class advection_diffusion_element : public element <datatype>
		{
		public:
			/*!*******************************************************************
			 * \param i_excess_0 The integer number of points evaluated in the adjacent element
			 * \param i_excess_n The integer number of points evaluated in the adjacent element
			 * \copydoc element::element ()
			 *********************************************************************/
			advection_diffusion_element (bases::axis *i_axis_n, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags);
			
			virtual ~advection_diffusion_element () {}
			
			virtual datatype calculate_timestep ();
		
		private:
			using element <datatype>::n;
			using element <datatype>::flags;
			using element <datatype>::name;
			using element <datatype>::normal_stream;
			using element <datatype>::cell;
			using element <datatype>::timestep;
			using element <datatype>::boundary_weights;
			using element <datatype>::inputParams;
			using element <datatype>::grids;
			using bases::element <datatype>::pointer;
			using element <datatype>::messenger_ptr;
		};
	} /* chebyshev */
	
	namespace fourier
	{
		/*!*******************************************************************
		 * \brief A Fourier implementation of the 1D element class
		 *********************************************************************/
		template <class datatype>
		class element : public one_d::element <datatype>
		{
		public:
			/*!*******************************************************************
			 * \copydoc one_d::element::element ()
			 *********************************************************************/
			element (bases::axis i_axis_n, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_axis_n, i_name, i_inputParams, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				one_d::element <datatype>::set_grid (new bases::fourier::grid <datatype> (axis_n, sqrt (2.0 / (n - 1.0)), messenger_ptr->linked (edge_0), messenger_ptr->linked (edge_n)));
				initialize (position);
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
			
			/*!*******************************************************************
			 * \copydoc one_d::element::initialize ()
			 *********************************************************************/
			virtual void initialize (int name, datatype* initial_conditions = NULL) {
				TRACE ("Initializing " << name);
				if (name == position && !initial_conditions) {
					one_d::element <datatype>::initialize (name, &(grids [0]->position ()));
				} else if (name == velocity && !initial_conditions){
					datatype scale = inputParams["init_cond_scale"].asDouble;
					datatype width = inputParams["init_cond_width"].asDouble;
					datatype mean = inputParams["init_cond_mean"].asDouble;
					datatype sigma = inputParams["init_cond_sigma"].asDouble;
					std::vector <datatype> init (n);
					datatype height, temp;
					height = std::max (scale * std::exp (- (width / 2.0 - mean) * (width / 2.0 - mean) / 2.0 / sigma / sigma), scale * std::exp (- (- width / 2.0 - mean) * (- width / 2.0 - mean) / 2.0 / sigma / sigma));
					for (int i = 0; i < n; ++i) {
						temp = scale * std::exp (- ((*this) (position, i) - mean) * ((*this) (position, i) - mean) / 2.0 / sigma / sigma) - height;
						if (temp > 0.0) {
							init [i] = temp;
						} else {
							init [i] = 0.0;
						}
					}
					one_d::element <datatype>::initialize (name, &init [0]);
				} else {
					one_d::element <datatype>::initialize (name, initial_conditions);
				}
				TRACE ("Initialized.");
			}
		
		protected:
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::axis_n;
			using one_d::element <datatype>::grids;
			using one_d::element <datatype>::inputParams;
			using one_d::element <datatype>::messenger_ptr;
		};
	} /* fourier */
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
