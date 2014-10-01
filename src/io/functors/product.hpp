/*!**********************************************************************
 * \file product.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PRODUCT_HPP_3123131F
#define PRODUCT_HPP_3123131F

#include "functor.hpp"

namespace io
{
	namespace functors
	{
		/*!**********************************************************************
		 * \brief Averages a two dimensional block of data in the horizontal direction
		 ************************************************************************/
		template <class datatype>
		class product_functor : public functor
		{
		private:
			datatype *data_1, *data_2; //!< A datatype pointer to the input data
			int n; //!< The integer horizontal extent of the data
			int m; //!< The integer vertical extent of the data
			std::vector <datatype> inner_data; //!< A vector of processed data to output

		public:
			/*!**********************************************************************
			 * \param i_data The datatype pointer to the data to average
			 * \param i_n The integer horizontal extent of the data
			 * \param i_m The integer vertical extent of the data
			 ************************************************************************/
			product_functor (int i_n, int i_m, datatype *i_data_1, datatype *i_data_2) : data_1 (i_data_1), data_2 (i_data_2), n (i_n), m (i_m) {
				inner_data.resize (n * m);
			}

			/*!**********************************************************************
			 * \brief Average the data and return a pointer to the first element
			 * 
			 * \return The first element of the averaged 1D array
			 ************************************************************************/
			void *calculate () {
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						inner_data [i * m + j] = data_1 [i * m + j] * data_2 [i * m + j];
					}
				}
				return &inner_data [0];
			}
		};
	} /* functors */
} /* io */

#endif /* end of include guard: PRODUCT_HPP_3123131F */
