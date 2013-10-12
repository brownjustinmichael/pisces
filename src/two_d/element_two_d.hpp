/*!**********************************************************************
 * \file element_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_TWO_D_HPP_CJ68F4IB
#define ELEMENT_TWO_D_HPP_CJ68F4IB

#include "../bases/element.hpp"

namespace two_d
{
	enum edges {
		edge_n0 = 0, // Start at 0, 0, increment by n
		edge_nn = 1, // Start at n, 0, increment by n
		edge_m0 = 2, // Start at 0, 0, increment by 1
		edge_mm = 3 // Start at 0, m, increment by 1
	}
	
	template <class datatype>
	class element
	{
	public:
		element (struct bases::axis i_axis_n, struct bases::axis i_axis_m, int i_name, io::parameter_map& i_inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 2, i_inputParams, i_messenger_ptr, i_flags),
		n (i_axis_n + 1),
		m (i_axis_m + 1) {
			positions [edge_n0] = i_axis_n.i_position_0;
			positions [edge_nn] = i_axis_n.i_position_n;
			positions [edge_m0] = i_axis_m.i_position_0;
			positions [edge_mm] = i_axis_m.i_position_n;
			
			excesses [edge_n0] = i_axis_n.i_excess_0;
			excesses [edge_nn] = i_axis_n.i_excess_n;
			excesses [edge_m0] = i_axis_m.i_excess_0;
			excesses [edge_mm] = i_axis_m.i_excess_n;
			
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
			
			cell.resize (n * m);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					cell_n [i + j * n] = i;
					cell_m [i + j * n] = j;
				}
			}
			
			std::ostringstream convert;
			convert << name;
			failsafe_dump.reset (new io::simple_output <datatype>  ("dump_" + convert.str () + ".dat", n));
			failsafe_dump->append (&cell_n [0]);
			failsafe_dump->append (&cell_m [0]);
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
		
		inline datatype& operator() (int name, int i, int j) {
			return bases::elelement <datatype>::operator() (name, i + j * n);
		}
		
		inline datatype* pointer (int name, int i, int j) {
			return bases::elelement <datatype>::pointer (name, i + j * n);
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual void initialize (int name, datatype* initial_conditions = NULL) {
			if (scalars [name].size () == (unsigned int) 0) {
				scalars [name].resize (n * m, 0.0);
			}
			if (initial_conditions) {
				utils::copy (n * m, initial_conditions, &(scalars [name]) [0]);
			}
			for (std::map <int, std::map <int, std::vector <datatype> > >::iterator iter = fixed_points.begin (); iter != fixed_points.end (); ++iter) {
				iter->second.resize (edge_size [iter->first]);
				utils::copy (edge_size [iter->first], pointer (name, edge_index [iter->first]), &(iter->second [name] [0]), edge_next [iter->first]);
			}
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element <datatype>::explicit_reset (); 
			for (iterator iter = bases::element <datatype>::begin (); iter != bases::element <datatype>::end (); ++iter) {
				if (*iter < 0) {
					utils::scale (n * m, 0.0, pointer (*iter));
				}
			}
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::execute_boundaries ()
		 *********************************************************************/
		inline void execute_boundaries () {
			for (std::map <int, std::map <int, std::vector <datatype> > >::iterator iter = fixed_points.begin (); iter != fixed_points.end (); ++iter) {
				for (std::map <int, std::vector <datatype> >::iterator i_name = iter->second.begin (); i_name != iter->second.end (); ++i_name) {
					utils::copy (edge_size [iter->first], &(i_name->second [0]), pointer (i_name->first, edge_index [iter->first]), 1, edge_next [iter->first]);
				}
			}
		}
	
	private:
		using bases::element <datatype>::name;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		
		int n; //!< The number of elements in each 1D array
		int m;
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
				element (struct bases::axis i_axis_n, struct bases::axis i_axis_m, int i_name, io::parameter_map& i_inputParams, bases::messenger* i_messenger_ptr, int i_flags) : 
				one_d::element <datatype> (i_axis_n, i_axis_m, i_name, i_inputParams, i_messenger_ptr, i_flags) {
					TRACE ("Instantiating...");
					initialize (position);
					two_d::element <datatype>::set_grid (new bases::fourier::grid <datatype> (n, n, sqrt (2.0 / (n - 1.0)), (*this) (position) - (*this) (position, n - 1)), 0);
					two_d::element <datatype>::set_grid (new bases::chebyshev::grid <datatype> (m, m, sqrt (2.0 / (m - 1.0)), (*this) (position) - (*this) (position, m - 1)), 1);
					TRACE ("Instantiated.");
				}
				virtual ~element () {}
			
				/*!*******************************************************************
				 * \copydoc one_d::element::initialize ()
				 *********************************************************************/
				virtual void initialize (int name, datatype* initial_conditions = NULL) {
					TRACE ("Initializing " << name);
					if (name == x_position && !initial_conditions) {
						std::vector <datatype> init (n);
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								init [i + j * n] = (i - excess_0) * (position_nn - position_n0) / (n - 1 - excess_nn - excess_n0) + position_0;
							}
						}
						two_d::element <datatype>::initialize (name, &init [0]);
					} else if (name == y_position && !initial_conditions) {
						std::vector <datatype> init (n);
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								init [i + j * n] = (j - excess_0) * (position_mm - position_m0) / (m - 1 - excess_mm - excess_m0) + position_0;
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
								init [i + j * n] = init [i];
							}
						}
						one_d::element <datatype>::initialize (name, &init [0]);
					} else {
						one_d::element <datatype>::initialize (name, initial_conditions);
					}
					TRACE ("Initialized.");
				}
		
			protected:
				using one_d::element <datatype>::position_n0;
				using one_d::element <datatype>::position_nn;
				using one_d::element <datatype>::position_m0;
				using one_d::element <datatype>::position_mm;
				using one_d::element <datatype>::excess_n0;
				using one_d::element <datatype>::excess_nn;
				using one_d::element <datatype>::excess_m0;
				using one_d::element <datatype>::excess_mm;
				using one_d::element <datatype>::n;
				using one_d::element <datatype>::m;
				using one_d::element <datatype>::inputParams;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
