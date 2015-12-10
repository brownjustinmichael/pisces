/*!**********************************************************************
 * \file slice.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SLICE_HPP_3123131F
#define SLICE_HPP_3123131F

#include <vector>
#include "functor.hpp"

namespace functors
{
	/*!**********************************************************************
	 * \brief Slices a two dimensional block of data horizontally
	 ************************************************************************/
	template <class datatype>
	class slice_functor : public functor
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		std::shared_ptr <functor> func; //!< A shared pointer to another functor object, for nesting
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		int i; //!< The integer location of the slice
		std::vector <datatype> inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_i The integer location of the slice
		 * \param i_data The datatype pointer to the data to average
		 ************************************************************************/
		slice_functor (int i_n, int i_m, int i_i, datatype *i_data) : data (i_data), n (i_n), m (i_m), i (i_i) {
			inner_data.resize (n);
		}
		
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_i The integer location of the slice
		 * \param i_func The shared pointer to the functor containing the processed data
		 ************************************************************************/
		slice_functor (int i_n, int i_m, int i_i, std::shared_ptr <functor> i_func) : data ((datatype *) i_func->calculate ()), func (i_func), n (i_n), m (i_m), i (i_i) {
			inner_data.resize (n);
		}

		/*!**********************************************************************
		 * \brief Average the data and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		void *calculate () {
			if (func) {
				func->calculate ();
			}
			for (int j = 0; j < n; ++j) {
				inner_data [j] = data [j * m + i];
			}
			return &inner_data [0];
		}
	};
} /* functors */

#endif /* end of include guard: SLICE_HPP_3123131F */
