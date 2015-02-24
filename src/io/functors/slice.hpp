/*!**********************************************************************
 * \file slice.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SLICE_HPP_3123131F
#define SLICE_HPP_3123131F

#include "functor.hpp"

namespace io
{
	namespace functors
	{
		/*!**********************************************************************
		 * \brief Averages a two dimensional block of data in the horizontal direction
		 ************************************************************************/
		template <class datatype>
		class slice_functor : public functor
		{
		private:
			datatype *data; //!< A datatype pointer to the input data
			std::shared_ptr <functor> func;
			int n; //!< The integer horizontal extent of the data
			int m; //!< The integer vertical extent of the data
			int i;
			std::vector <datatype> inner_data; //!< A vector of processed data to output

		public:
			/*!**********************************************************************
			 * \param i_data The datatype pointer to the data to average
			 * \param i_n The integer horizontal extent of the data
			 * \param i_m The integer vertical extent of the data
			 ************************************************************************/
			slice_functor (int i_n, int i_m, int i_i, datatype *i_data) : data (i_data), n (i_n), m (i_m), i (i_i) {
				inner_data.resize (n);
			}
			
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
} /* io */

#endif /* end of include guard: SLICE_HPP_3123131F */
