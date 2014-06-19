/*!***********************************************************************
 * \file element_one_d.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_3SURDTOH
#define ELEMENT_HPP_3SURDTOH

#include "../utils/messenger.hpp"
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include "transform_one_d.hpp"
#include "../bases/element.hpp"
#include "../utils/utils.hpp"
#include "../utils/rezone.hpp"
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
	 * storage and indexing facilities. The plans should be added in a 
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
		element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 1, i_params, i_messenger_ptr, i_flags),
		n (i_axis_n.get_n ()) {
			axes [0] = i_axis_n;
			
			if (i_messenger_ptr->get_id () == 0) {
				alpha_0 = 0.0;
			} else {
				alpha_0 = 0.5;
			}
			if (i_messenger_ptr->get_id () == i_messenger_ptr->get_np () - 1) {
				alpha_n = 0.0;
			} else {
				alpha_n = 0.5;
			}
			
			max_timestep = i_params.get <datatype> ("time.max");
		}
		
		virtual ~element () {
			// printf ("Destroying one_d element\n");
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			
			scalars [name].resize ((*grids [0]).get_ld (), 0.0);
			if (name == position) {
				utils::copy (n, &((*grids [0]) [0]), this->ptr (name));
			} else if (initial_conditions) {
				utils::copy (n, initial_conditions, this->ptr (name));
			}
			
			TRACE ("Initialized " << name << ".");
			return this->ptr (name);
		}
		
		virtual datatype calculate_timestep (int i, io::virtual_dump *dump = NULL) = 0;

		inline datatype calculate_min_timestep (io::virtual_dump *dump = NULL) {
			datatype shared_min = max_timestep;
			#pragma omp parallel 
			{
				datatype min = std::numeric_limits <datatype>::max ();
				#pragma omp for nowait
					for (int i = 1; i < (dump ? dump->dims ["z"] [1] : n) - 1; ++i) {
						min = std::min (calculate_timestep (i, dump), min);
					}
				#pragma omp critical 
				{
					shared_min = std::min (shared_min, min);
				}
			}
			if (shared_min < timestep || shared_min > 2.0 * timestep) {
				return shared_min * 0.8;
			} else {
				return timestep;
			}
		}
		
		virtual io::virtual_dump *make_dump (int flags = 0x00) {
			std::shared_ptr <io::output> virtual_output (new io::formatted_output <io::formats::two_d::virtual_format> ("one_d/dump", n));
			bases::element <datatype>::setup_output (virtual_output);
			
			virtual_output->to_file ();
			return &io::virtual_dumps ["one_d/dump"];
		}
		
		virtual std::shared_ptr <bases::grid <datatype>> generate_grid (bases::axis *axis, int index = -1) = 0;
		
		virtual io::virtual_dump *make_rezoned_dump (datatype *positions, io::virtual_dump *old_dump, int flags = 0x00) {
			bases::axis vertical_axis (n, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () ? 0 : 1);
			std::shared_ptr <bases::grid <datatype>> vertical_grid = generate_grid (&vertical_axis);
			
			utils::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, old_dump, &io::virtual_dumps ["one_d/new_dump"]);
			
			return &io::virtual_dumps ["one_d/new_dump"];
		}
		
		/*
			TODO Combine these?
		*/
		
		virtual void get_zoning_positions (datatype *positions) {
			if (messenger_ptr->get_id () == 0) {
				datatype temp [messenger_ptr->get_np () * 2];
				temp [0] = axes [0].get_position_0 ();
				temp [1] = axes [0].get_position_n ();
				messenger_ptr->template gather <datatype> (2, temp, temp);
				for (int i = 0; i < messenger_ptr->get_np (); ++i) {
					positions [i] = temp [2 * i];
				}
				positions [messenger_ptr->get_np ()] = temp [messenger_ptr->get_np () * 2 - 1];
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			} else {
				datatype temp [2] = {axes [0].get_position_0 (), axes [0].get_position_n ()};
				messenger_ptr->template gather <datatype> (2, temp, temp);
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			}
		}
		
	protected:
		using bases::element <datatype>::scalars;
		using bases::element <datatype>::params;
		using bases::element <datatype>::grids;
		using bases::element <datatype>::name;
		using bases::element <datatype>::messenger_ptr;
		using bases::element <datatype>::ptr;
		using bases::element <datatype>::timestep;
		using bases::element <datatype>::axes;
		
		bases::axis *axis_n;
		datatype alpha_0, alpha_n, max_timestep;
		int n;
		std::vector<int> cell; //!< An integer array for tracking each cell number for output
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
			element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_axis_n, i_name, i_params, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				// one_d::element <datatype>::set_grid (new bases::chebyshev::grid <datatype> (axis_n));
				initialize (position, "x");
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
			
			const int &get_mode () {
				return mode;
			}
			
			virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
				TRACE ("Initializing " << name << "...");
				one_d::element <datatype>::_initialize (name, initial_conditions, flags);
				if (!(flags & no_transform) && (name != position)) {
#ifdef _CUDA
					element <datatype>::add_transform (name, std::shared_ptr <cuda::master_transform <datatype> > (new cuda::master_transform <datatype> (*grids [0], ptr (name), NULL, forward_vertical | inverse_vertical, &element_flags [state], &element_flags [name])));
#else
					element <datatype>::add_transform (name, std::shared_ptr <master_transform <datatype> > (new master_transform <datatype> (*grids [0], ptr (name), NULL, forward_vertical | inverse_vertical, &element_flags [state], &element_flags [name])));
#endif /* _CUDA */
				}
				return this->ptr (name);
			}
			
			virtual std::shared_ptr <bases::grid <datatype>> generate_grid (bases::axis *axis, int index = -1) {
				if (index == 0 || index == -1) {
					return std::shared_ptr <bases::grid <datatype>> (new bases::chebyshev::grid <datatype> (axis));
				} else {
					throw 0;
				}
			}
			
			using one_d::element <datatype>::initialize;
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::grids;
			using one_d::element <datatype>::params;
			using one_d::element <datatype>::messenger_ptr;
			using one_d::element <datatype>::ptr;
			using one_d::element <datatype>::element_flags;
			
			static int mode;
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
			nonlinear_diffusion_element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_flags);
		
			virtual ~nonlinear_diffusion_element () {}
		
			virtual datatype calculate_timestep (int i, io::virtual_dump *dump);
	
			using element <datatype>::initialize;
			using element <datatype>::n;
			using element <datatype>::element_flags;
			using element <datatype>::name;
			using element <datatype>::cell;
			using element <datatype>::timestep;
			using element <datatype>::params;
			using element <datatype>::grids;
			using bases::element <datatype>::ptr;
			using bases::element <datatype>::matrix_ptr;
			using element <datatype>::messenger_ptr;
			using element <datatype>::alpha_0;
			using element <datatype>::alpha_n;
			using element <datatype>::solvers;
			
		private:
			datatype advection_coeff, cfl;
			datatype *position_ptr, *velocity_ptr;
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
			advection_diffusion_element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_flags);
		
			virtual ~advection_diffusion_element () {}
		
			virtual datatype calculate_timestep (int i, io::virtual_dump *dump);
	
			using element <datatype>::initialize;
			using element <datatype>::n;
			using element <datatype>::element_flags;
			using element <datatype>::name;
			using element <datatype>::cell;
			using element <datatype>::timestep;
			using bases::element <datatype>::params;
			using element <datatype>::grids;
			using bases::element <datatype>::ptr;
			using bases::element <datatype>::matrix_ptr;
			using element <datatype>::messenger_ptr;
			using element <datatype>::alpha_0;
			using element <datatype>::alpha_n;
			using element <datatype>::solvers;
		
		private:
			datatype advection_coeff, cfl;
			datatype *position_ptr, *velocity_ptr;
		};
	} /* chebyshev */
	
	namespace cosine
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
			element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_axis_n, i_name, i_params, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				initialize (position, "x");
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
			
			const int &get_mode () {
				return mode;
			}
		
			virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
				TRACE ("Initializing " << name << "...");
				one_d::element <datatype>::_initialize (name, initial_conditions, flags);
				if (!(flags & no_transform) && (name != position)) {
#ifdef _CUDA
					element <datatype>::add_transform (name, std::shared_ptr <cuda::master_transform <datatype> > (new cuda::master_transform <datatype> (*grids [0], ptr (name), NULL, forward_vertical | inverse_vertical, &element_flags [state], &element_flags [name])));
#else
					element <datatype>::add_transform (name, std::shared_ptr <master_transform <datatype> > (new master_transform <datatype> (*grids [0], ptr (name), NULL, forward_vertical | inverse_vertical, &element_flags [state], &element_flags [name])));
#endif /* _CUDA */
					
				}
				return this->ptr (name);
			}
			
			virtual std::shared_ptr <bases::grid <datatype>> generate_grid (bases::axis *axis, int index = -1) {
				if (index == 0 || index == -1) {
					return std::shared_ptr <bases::grid <datatype>> (new bases::cosine::grid <datatype> (axis));
				} else {
					throw 0;
				}
			}
		
			using one_d::element <datatype>::initialize;
			using one_d::element <datatype>::element_flags;
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::grids;
			using one_d::element <datatype>::params;
			using one_d::element <datatype>::messenger_ptr;
			using one_d::element <datatype>::solvers;
			using one_d::element <datatype>::ptr;
			
			static int mode;
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
			advection_diffusion_element (bases::axis i_axis_n, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_flags);
	
			virtual ~advection_diffusion_element () {}
	
			virtual datatype calculate_timestep (int i, io::virtual_dump *dump);

			using element <datatype>::initialize;
			using bases::element <datatype>::ptr;
			using bases::element <datatype>::matrix_ptr;
			
		private:
			using element <datatype>::n;
			using element <datatype>::element_flags;
			using element <datatype>::name;
			using element <datatype>::cell;
			using element <datatype>::grids;
			using element <datatype>::messenger_ptr;
			using element <datatype>::timestep;
			using bases::element <datatype>::params;
			using element <datatype>::alpha_0;
			using element <datatype>::alpha_n;
			using element <datatype>::solvers;
			datatype advection_coeff, cfl;
			datatype *position_ptr, *velocity_ptr;
		};
		/*
			TODO Fix member permissions
		*/
	} /* cosine */
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
