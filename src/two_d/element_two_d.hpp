/*!**********************************************************************
 * \file element_two_d.hpp
 * /Users/justinbrown/Dropbox/spectral_element/src
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
		edge_00 = 0, // Start at 0, 0, increment by n
		edge_nm = 1, // Start at 0, m, increment by n
		edge_0m = 2, // Start at 0, 0, increment by 1
		edge_n0 = 3 // Start at n, 0, increment by 1
	}
	
	template <class datatype>
	class element
	{
	public:
		element (int i_n, datatype i_position_00, datatype i_position_nm, int i_m, datatype i_position_0m, datatype i_position_n0, int i_name, io::parameter_map& i_inputParams, bases::messenger <datatype>* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 2, i_inputParams, i_messenger_ptr, i_flags),
		n (i_n + 1),
		m (i_m + 1) {
			positions.resize (4);
			positions [edge_00] = i_position_00;
			positions [edge_nm] = i_position_nm;
			positions [edge_0m] = i_position_0m;
			positions [edge_n0] = i_position_n0;
			
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
				bases::element <datatype>::add_name (name);
				scalars [name].resize (n * m, 0.0);
			}
			if (initial_conditions) {
				utils::copy (n * m, initial_conditions, &(scalars [name]) [0]);
			}
			fixed_points_00 [name].resize (m);
			fixed_points_nm [name].resize (m);
			fixed_points_0m [name].resize (n);
			fixed_points_n0 [name].resize (n);
			
			utils::copy (m, pointer (name, 0, 0), &(fixed_points_00 [name] [0]), n);
			utils::copy (m, pointer (name, 0, m), &(fixed_points_nm [name] [0]), n);
			utils::copy (n, pointer (name, 0, 0), &(fixed_points_0m [name] [0]));
			utils::copy (n, pointer (name, n, 0), &(fixed_points_n0 [name] [0]));
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
			for (iterator iter = bases::element <datatype>::begin (); iter != bases::element <datatype>::end (); ++iter) {
				if (!(messenger_ptr->linked (edge_00))) {
					utils::copy (m, &(fixed_points_00 [name] [0]), pointer (name, 0, 0), 1, n);
				}
				if (!(messenger_ptr->linked (edge_nm))) {
					utils::copy (m, &(fixed_points_nm [name] [0]), pointer (name, 0, m), 1, n);
				}
				if (!(messenger_ptr->linked (edge_0m))) {
					utils::copy (m, &(fixed_points_0m [name] [0]), pointer (name, 0, 0));
				}
				if (!(messenger_ptr->linked (edge_n0))) {
					utils::copy (m, &(fixed_points_n0 [name] [0]), pointer (name, n, 0));
				}
			}
		}
	
	private:
		using bases::element <datatype>::name;
		using bases::element <datatype>::names;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		typedef typename bases::element <datatype>::iterator iterator;
		
		int n; //!< The number of elements in each 1D array
		int m;
		std::vector <datatype> positions; //!< A vector of the edge positions
		std::vector<int> cell_n; //!< An integer array for tracking each cell number for output
		std::vector<int> cell_m; //!< An integer array for tracking each cell number for output

		std::map <int, std::vector <datatype> > scalars; //!< A vector of scalar vectors
		
		std::map <int, int [2]> edge_map;
		std::map <int, std::map <int, std::vector <datatype> > > fixed_points;
		std::map <int, std::vector <datatype> > fixed_points_00; //!< The initial values of the scalars at index 0
		std::map <int, std::vector <datatype> > fixed_points_nm; //!< The initial values of the scalars at index 0
		std::map <int, std::vector <datatype> > fixed_points_n0; //!< The initial values of the scalars at index 0
		std::map <int, std::vector <datatype> > fixed_points_0m; //!< The initial values of the scalars at index 0
	};
	
	namespace fourier
	{
		namespace chebyshev
		{
			
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
