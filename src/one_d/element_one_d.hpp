/*!***********************************************************************
 * \file element_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_3SURDTOH
#define ELEMENT_HPP_3SURDTOH

#include "../bases/messenger.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include "transform_one_d.hpp"
#include "../bases/element.hpp"
#include "../utils/utils.hpp"
#include "../bases/grid.hpp"
#include "../config.hpp"
	
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
	
	enum initialize_flags {
		uniform_n = 0x01
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
		element (bases::axis *i_axis_n, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 1, i_params, i_messenger_ptr, i_flags),
		axis_n (i_axis_n),
		n (axis_n->n) {
			cell.resize (n);
			for (int i = 0; i < n; ++i) {
				cell [i] = i;
			}
			
			edge_map [edge_0] = 0;
			edge_map [edge_n] = n - 1;
			
			if (i_messenger_ptr->get_id () == 0) {
				fixed_points [edge_0];
				alpha_0 = 0.0;
			} else {
				alpha_0 = 0.5;
			}
			if (i_messenger_ptr->get_id () == i_messenger_ptr->get_np () - 1) {
				fixed_points [edge_n];
				alpha_n = 0.0;
			} else {
				alpha_n = 0.5;
			}
			
			std::ostringstream convert;
			convert << name;
			failsafe_dump.reset (new io::output (new io::one_d::ascii (n), "dump.dat"));
			failsafe_dump->append ("i", &cell [0]);
		}
		
		virtual ~element () {
			// printf ("Destroying one_d element\n");
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual void initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			
			scalars [name].resize (n, 0.0);
			if (name == position) {
				initial_conditions = &(grids [0]->position ());
			}
			if (initial_conditions) {
				utils::copy (n, initial_conditions, this->ptr (name));
			}
			for (std::map <int, int>::iterator i_edge = edge_map.begin (); i_edge != edge_map.end (); ++i_edge) {
				fixed_points [i_edge->first] [name] = scalars [name] [i_edge->second];
			}
			std::stringstream stream;
			stream << name;
			failsafe_dump->template append <datatype> (stream.str (), &(scalars [name]) [0]);

			TRACE ("Initialized " << name << ".");
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
		
	protected:
		using bases::element <datatype>::scalars;
		using bases::element <datatype>::grids;
		using bases::element <datatype>::name;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		using bases::element <datatype>::ptr;
		
		bases::axis *axis_n;
		datatype alpha_0, alpha_n;
		int &n;
		std::vector<int> cell; //!< An integer array for tracking each cell number for output
		
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
			element (bases::axis *i_axis_n, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_axis_n, i_name, i_params, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				one_d::element <datatype>::set_grid (new bases::chebyshev::grid <datatype> (axis_n));
				initialize (position);
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
			
			virtual void initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
				TRACE ("Initializing " << name << "...");
				one_d::element <datatype>::initialize (name, initial_conditions, flags);
				if (!(flags & no_transform) && (name != position)) {
					element <datatype>::add_inverse_vertical_transform (new fftw_cosine <datatype> (*grids [0], ptr (name)));
					element <datatype>::add_forward_vertical_transform (new fftw_cosine <datatype> (*grids [0], ptr (name)));
				}
			}
			
		protected:
			using one_d::element <datatype>::initialize;
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::axis_n;
			using one_d::element <datatype>::grids;
			using one_d::element <datatype>::params;
			using one_d::element <datatype>::messenger_ptr;
			using one_d::element <datatype>::ptr;
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
			advection_diffusion_element (bases::axis *i_axis_n, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags);
			
			virtual ~advection_diffusion_element () {}
			
			virtual datatype calculate_timestep ();
		
		private:
			using element <datatype>::initialize;
			using element <datatype>::n;
			using element <datatype>::element_flags;
			using element <datatype>::name;
			using element <datatype>::normal_stream;
			using element <datatype>::cell;
			using element <datatype>::timestep;
			using element <datatype>::params;
			using element <datatype>::grids;
			using bases::element <datatype>::ptr;
			using bases::element <datatype>::matrix_ptr;
			using element <datatype>::messenger_ptr;
			using element <datatype>::alpha_0;
			using element <datatype>::alpha_n;
		};
		
		/*!*******************************************************************
		 * \brief A simple implementation of the element class with diffusion
		 * 
		 * This class contains a full element's capacity to run a single 
		 * element diffusion in 1D with constant timestep.
		 *********************************************************************/
		template <class datatype>
		class nonlinear_diffusion_element : public element <datatype>
		{
		public:
			/*!*******************************************************************
			 * \param i_excess_0 The integer number of points evaluated in the adjacent element
			 * \param i_excess_n The integer number of points evaluated in the adjacent element
			 * \copydoc element::element ()
			 *********************************************************************/
			nonlinear_diffusion_element (bases::axis *i_axis_n, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags);
		
			virtual ~nonlinear_diffusion_element () {}
		
			virtual datatype calculate_timestep ();
	
		private:
			using element <datatype>::initialize;
			using element <datatype>::n;
			using element <datatype>::element_flags;
			using element <datatype>::name;
			using element <datatype>::normal_stream;
			using element <datatype>::cell;
			using element <datatype>::timestep;
			using element <datatype>::params;
			using element <datatype>::grids;
			using bases::element <datatype>::ptr;
			using bases::element <datatype>::matrix_ptr;
			using element <datatype>::messenger_ptr;
			using element <datatype>::alpha_0;
			using element <datatype>::alpha_n;
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
			element (bases::axis i_axis_n, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_axis_n, i_name, i_params, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				one_d::element <datatype>::set_grid (new bases::fourier::grid <datatype> (axis_n));
				initialize (position);
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
		
		protected:
			using one_d::element <datatype>::initialize;
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::axis_n;
			using one_d::element <datatype>::grids;
			using one_d::element <datatype>::params;
			using one_d::element <datatype>::messenger_ptr;
		};
	} /* fourier */
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
