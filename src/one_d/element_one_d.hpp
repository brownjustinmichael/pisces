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
#include <vector>
#include <map>
#include "../config.hpp"
#include "../bases/element.hpp"
#include "../bases/plan.hpp"
#include "../utils/utils.hpp"
#include "../utils/chebyshev.hpp"
	
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
		element (int i_n, datatype i_position_0, datatype i_position_n, int i_name, io::parameter_map& i_inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 1, i_inputParams, i_messenger_ptr, i_flags) {
			n = i_n;
			position_0 = i_position_0;
			position_n = i_position_n;

			cell.resize (i_n);
			for (int i = 0; i < i_n; ++i) {
				cell [i] = i;
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
				bases::element <datatype>::add_name (name);
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
		 * \copydoc bases::element <datatype>::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element <datatype>::explicit_reset (); 
			for (iterator iter = bases::element <datatype>::begin (); iter != bases::element <datatype>::end (); ++iter) {
				if (*iter < 0) {
					utils::scale (n, 0.0, &((*this) [*iter]));
				}
			}
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::execute_boundaries ()
		 *********************************************************************/
		inline void execute_boundaries () {
			if (messenger_ptr->linked (edge_0)) {
				for (iterator iter = bases::element <datatype>::begin (); iter != bases::element <datatype>::end (); ++iter) {
					(*this) (*iter, 0) = fixed_points_0 [*iter];
				}
			}
			if (messenger_ptr->linked (edge_n)) {
				for (iterator iter = bases::element <datatype>::begin (); iter != bases::element <datatype>::end (); ++iter) {
					(*this) (*iter, n - 1) = fixed_points_n [*iter];
				}
			}
			
			/*
				TODO Not sure why these aren't !linked...
			*/
		}
		
	protected:
		using bases::element <datatype>::name;
		using bases::element <datatype>::names;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		typedef typename bases::element <datatype>::iterator iterator;
		
		int n; //!< The number of elements in each 1D array
		datatype position_0; //!< The datatype position of index 0
		datatype position_n; //!< The datatype position of index n - 1
		std::vector<int> cell; //!< An integer array for tracking each cell number for output

		std::map <int, std::vector <datatype>> scalars; //!< A vector of scalar vectors
		std::map <int, datatype> fixed_points_0; //!< The initial values of the scalars at index 0
		std::map <int, datatype> fixed_points_n; //!< The initial values of the scalars at index n - 1
		
		/*
			TODO Perhaps there's a better way to handle the fixed points
		*/
	};

	namespace chebyshev
	{
		/*!**********************************************************************
		 * \brief Helper to calculate positions in chebyshev elements
		 * 
		 * \param n The integer number of data points
		 * \param i The integer element index for the calculation
		 * \param excess_0 The integer number of excess data points on the edge_0 side
		 * \param position_0 The datatype position at index excess_0
		 * \param excess_n The integer number of excess data points on the edge_n side
		 * \param position_n The datatype position at index n - 1 - excess_n
		 * 
		 * \return The datatype position of the given index
		 ************************************************************************/
		template <class datatype>
		datatype return_position (int n, int i, int excess_0, datatype position_0, int excess_n, datatype position_n) {
			datatype pioN = std::acos (-1.0) / (n - 1);
			datatype scale = (position_0 - position_n) / (std::cos (excess_0 * pioN) - std::cos ((n - 1 - excess_n) * pioN));
			datatype initial = position_0 - scale * std::cos (excess_0 * pioN);
			return scale * std::cos (i * pioN) + initial;
		}
		
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
			element (int i_n, datatype i_position_0, datatype i_position_n, int i_name, io::parameter_map& i_inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags) : 
			one_d::element <datatype> (i_n, i_position_0, i_position_n, i_name, i_inputParams, i_messenger_ptr, i_flags) {
				TRACE ("Instantiating...");
				initialize (position);
				one_d::element <datatype>::set_grid (new chebyshev_grid <datatype> (i_n, i_n, sqrt (2.0 / (i_n - 1.0)), position_0 - position_n));
				TRACE ("Instantiated.");
			}
			virtual ~element () {}
				
			/*!*******************************************************************
			 * \copydoc one_d::element::initialize ()
			 *********************************************************************/
			virtual void initialize (int name, datatype* initial_conditions = NULL) {
				TRACE ("Initializing " << name);
				if (name == position && !initial_conditions) {
					datatype pioN = std::acos (-1.0) / (n - 1);
					datatype scale = (position_0 - position_n) / 2.0;
					datatype initial_position = (position_0 + position_n) / 2.0;
					std::vector <datatype> init (n);
					for (int i = 0; i < n; ++i) {
						init [i] = scale * std::cos (i * pioN) + initial_position;
					}
					one_d::element <datatype>::initialize (name, &init [0]);
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
			using one_d::element <datatype>::position_0;
			using one_d::element <datatype>::position_n;
			using one_d::element <datatype>::n;
			using one_d::element <datatype>::inputParams;
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
			advection_diffusion_element (int i_n, datatype i_position_0, datatype i_position_n, int i_excess_0, int i_excess_n, int i_name, io::parameter_map& i_inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags);
			
			virtual ~advection_diffusion_element () {}
		
			/*!*******************************************************************
			 * \copydoc element::implicit_reset ()
			 *********************************************************************/
			inline void implicit_reset () {
				element <datatype>::implicit_reset ();
			
				if (!(flags & factorized)) {
					utils::scale (n * n, 0.0, &matrix [0]);
				}
			}
			
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
			using element <datatype>::grid;
			typedef typename element <datatype>::iterator iterator;
		
			std::vector<datatype> matrix; //!< A vector containing the datatype matrix used in the implicit solver
			std::vector<datatype> temp_matrix; //!< A vector containing the datatype matrix used in the implicit solver
		};
	} /* chebyshev */
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
