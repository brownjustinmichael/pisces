/*!**********************************************************************
 * \file element_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_TWO_D_HPP_CJ68F4IB
#define ELEMENT_TWO_D_HPP_CJ68F4IB

#include "../config.hpp"
#include "../bases/element.hpp"
#include <cmath>

namespace two_d
{
	enum edges {
		edge_n0 = 0, // Start at 0, 0, increment by n
		edge_nn = 1, // Start at n, 0, increment by n
		edge_m0 = 2, // Start at 0, 0, increment by 1
		edge_mm = 3 // Start at 0, m, increment by 1
	};
	
	template <class datatype>
	class element : public bases::element <datatype>
	{
	public:
		element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 2, i_inputParams, i_messenger_ptr, i_flags),
		axis_n (i_axis_n),
		axis_m (i_axis_m),
		n (axis_n->n), m (axis_m->n) {
			TRACE ("Instantiating...");
			
			positions [edge_n0] = i_axis_n->position_0;
			positions [edge_nn] = i_axis_n->position_n;
			positions [edge_m0] = i_axis_m->position_0;
			positions [edge_mm] = i_axis_m->position_n;
			
			excesses [edge_n0] = i_axis_n->excess_0;
			excesses [edge_nn] = i_axis_n->excess_n;
			excesses [edge_m0] = i_axis_m->excess_0;
			excesses [edge_mm] = i_axis_m->excess_n;
			
			edge_index [edge_n0] = 0;
			edge_next [edge_n0] = n;
			edge_size [edge_n0] = m;
			
			edge_index [edge_nn] = n - 1;
			edge_next [edge_nn] = n;
			edge_size [edge_nn] = m;
			
			edge_index [edge_m0] = 0;
			edge_next [edge_m0] = 1;
			edge_size [edge_m0] = n;
			
			edge_index [edge_mm] = (m - 1) * n;
			edge_next [edge_mm] = 1;
			edge_size [edge_mm] = n;
			
			cell_n.resize (n);
			cell_m.resize (m);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					cell_n [i * m + j] = i;
					cell_m [i * m + j] = j;
				}
			}
			
			std::ostringstream convert;
			convert << name;
			failsafe_dump.reset (new io::simple_output <datatype>  ("dump_" + convert.str () + ".dat", n));
			failsafe_dump->append (&cell_n [0]);
			failsafe_dump->append (&cell_m [0]);
			
			TRACE ("Instantiated.");
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
		
		inline datatype& operator() (int name, int i = 0, int j = 0) {
			return bases::element <datatype>::operator() (name, i * m + j);
		}
		
		inline datatype* pointer (int name, int i = 0, int j = 0) {
			return bases::element <datatype>::pointer (name, i * m + j);
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual void initialize (int name, datatype* initial_conditions = NULL) {
			TRACE ("Initializing...");
			
			if (scalars [name].size () == (unsigned int) 0) {
				scalars [name].resize (n * m, 0.0);
			}
			if (initial_conditions) {
				utils::copy (n * m, initial_conditions, &(scalars [name]) [0]);
			}
			for (typename std::map <int, std::map <int, std::vector <datatype> > >::iterator iter = fixed_points.begin (); iter != fixed_points.end (); ++iter) {
				iter->second [name].resize (edge_size [iter->first]);
				utils::copy (edge_size [iter->first], pointer (name, edge_index [iter->first]), &(iter->second [name] [0]), edge_next [iter->first]);
			}
			
			TRACE ("Initialized.");
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element <datatype>::explicit_reset ();
			
			TRACE ("Explicit reset start...");
			
			for (typename std::map <int, std::vector <datatype> >::iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				if (iter->first < 0) {
					utils::scale (n * m, 0.0, &(iter->second) [0]);
				}
			}
			
			TRACE ("Explicit reset end.");
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::execute_boundaries ()
		 *********************************************************************/
		inline void execute_boundaries () {
			TRACE ("Executing boundaries...");
			
			for (typename std::map <int, std::map <int, std::vector <datatype> > >::iterator iter = fixed_points.begin (); iter != fixed_points.end (); ++iter) {
				for (typename std::map <int, std::vector <datatype> >::iterator i_name = iter->second.begin (); i_name != iter->second.end (); ++i_name) {
					utils::copy (edge_size [iter->first], &(i_name->second [0]), pointer (i_name->first, edge_index [iter->first]), 1, edge_next [iter->first]);
				}
			}
			
			TRACE ("Boundaries executed.");
		}
	
	protected:
		using bases::element <datatype>::name;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		
		bases::axis *axis_n, *axis_m;
		int &n; //!< The number of elements in each 1D array
		int &m;
		std::map <int, datatype> positions; //!< A vector of the edge positions
		std::map <int, int> excesses; //!< A vector of the edge positions
		std::vector<int> cell_n; //!< An integer array for tracking each cell number for output
		std::vector<int> cell_m; //!< An integer array for tracking each cell number for output

		std::map <int, std::vector <datatype> > scalars; //!< A vector of scalar vectors
		
		std::map <int, int> edge_index;
		std::map <int, int> edge_next;
		std::map <int, int> edge_size;
		std::map <int, std::map <int, std::vector <datatype> > > fixed_points;
	};
	
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class element : public two_d::element <datatype>
			{
			public:
				/*!*******************************************************************
				 * \copydoc one_d::element::element ()
				 *********************************************************************/
				element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
				two_d::element <datatype> (i_axis_n, i_axis_m, i_name, i_inputParams, i_messenger_ptr, i_flags) {
					TRACE ("Instantiating...");
					
					two_d::element <datatype>::set_grid (new bases::fourier::grid <datatype> (axis_n, sqrt (2.0 / (n - 1.0))), 0);
					two_d::element <datatype>::set_grid (new bases::chebyshev::grid <datatype> (axis_m, sqrt (2.0 / (m - 1.0))), 1);
					initialize (x_position);
					initialize (z_position);
					
					TRACE ("Instantiated.");
				}
				virtual ~element () {}
			
				/*!*******************************************************************
				 * \copydoc one_d::element::initialize ()
				 *********************************************************************/
				virtual void initialize (int name, datatype* initial_conditions = NULL) {
					
					TRACE ("Initializing " << name << "...");
					
					if (name == x_position && !initial_conditions) {
						std::vector <datatype> init (n * m);
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								init [i * m + j] = (i - excesses [edge_n0]) * (positions [edge_nn] - positions [edge_n0]) / (n - 1 - excesses [edge_nn] - excesses [edge_n0]) + positions [edge_n0];
							}
						}
						two_d::element <datatype>::initialize (name, &init [0]);
					} else if (name == z_position && !initial_conditions) {
						std::vector <datatype> init (n * m);
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								init [i * m + j] = (j - excesses [edge_m0]) * (positions [edge_mm] - positions [edge_m0]) / (m - 1 - excesses [edge_mm] - excesses [edge_m0]) + positions [edge_m0];
							}
						}
						two_d::element <datatype>::initialize (name, &init [0]);
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
							for (int j = 0; j < m; ++j) {
								init [i * m + j] = init [i];
							}
						}
						two_d::element <datatype>::initialize (name, &init [0]);
					} else {
						two_d::element <datatype>::initialize (name, initial_conditions);
					}
					
					TRACE ("Initialized.");
				}
		
			protected:
				using two_d::element <datatype>::positions;
				using two_d::element <datatype>::excesses;
				using two_d::element <datatype>::n;
				using two_d::element <datatype>::m;
				using two_d::element <datatype>::axis_n;
				using two_d::element <datatype>::axis_m;
				using two_d::element <datatype>::inputParams;
			};
			
			class advection_diffusion_element
			{
			public:
				advection_diffusion_element (struct bases::axis i_axis_n, struct bases::axis i_axis_m, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags);
				
				virtual ~advection_diffusion_element () {}
			
			private:
				/* data */
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
