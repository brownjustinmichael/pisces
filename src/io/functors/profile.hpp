/*!**********************************************************************
 * \file profile.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PROFILE_HPP_3123131F
#define PROFILE_HPP_3123131F

#include <vector>
#include "functor.hpp"

namespace functors
{
	/*!**********************************************************************
	 * \brief Profiles a two dimensional block of data vertically
	 ************************************************************************/
	template <class datatype>
	class profile_functor : public functor
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		std::shared_ptr <functor> func; //!< A shared pointer to another functor object, for nesting
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_i The integer location of the profile
		 * \param i_data The datatype pointer to the data to average
		 ************************************************************************/
		profile_functor (int i_n, int i_m, datatype *i_data) : data (i_data), n (i_n), m (i_m) {
			inner_data.resize (m);
		}
		
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_i The integer location of the profile
		 * \param i_func The shared pointer to the functor containing the processed data
		 ************************************************************************/
		profile_functor (int i_n, int i_m, std::shared_ptr <functor> i_func) : data ((datatype *) i_func->calculate ()), func (i_func), n (i_n), m (i_m) {
			DEBUG(" " << data);
			DEBUG(" " << func);
			inner_data.resize (m);
		}

		/*!**********************************************************************
		 * \brief Average the data and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		void *calculate () {
			DEBUG("Cacluating");
			if (func) {
				DEBUG("Cacluating func");
				func->calculate ();
			}
			for (int j = 0; j < m; ++j) {
				inner_data[j] = 0.0;
				for (int i = 0; i < n; ++i) {
					inner_data [j] += data [i * m + j];
				}
				inner_data[j] /= n;
			}
			return &inner_data [0];
		}
	};
} /* functors */

#endif /* end of include guard: PROFILE_HPP_3123131F */
