/*!**********************************************************************
 * \file div.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIV_HPP_A5256DE8
#define DIV_HPP_A5256DE8

#include <vector>
#include <cmath>

#include "functor.hpp"

namespace functors
{
	/*!**********************************************************************
	 * \brief Takes the divergence of a two dimensional vector field
	 ************************************************************************/
	template <class datatype>
	class div_functor : public functor
	{
	private:
		datatype *pos_x; //!< A datatype pointer to the horizontal position data
		datatype *pos_z; //!< A datatype pointer to the vertical position data
		datatype *data_x; //!< A datatype pointer to the x component of the vector data
		datatype *data_z; //!< A datatype pointer to the z component of the vector data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_pos_x A datatype pointer to the horizontal position data
		 * \param i_pos_z A datatype pointer to the vertical position data
		 * \param i_data_x A datatype pointer to the x component of the vector data
		 * \param i_data_z A datatype pointer to the z component of the vector data
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		div_functor (datatype *i_pos_x, datatype *i_pos_z, datatype *i_data_x, datatype *i_data_z, int i_n, int i_m) : pos_x (i_pos_x), pos_z (i_pos_z), data_x (i_data_x), data_z (i_data_z), n (i_n), m (i_m) {
			inner_data.resize (n * m);
		}

		/*!**********************************************************************
		 * \brief Take the divergence of the data and return a pointer to the first element
		 * 
		 * \return The first element of the 2D divergence array
		 ************************************************************************/
		void *calculate () {
			// Calculate the vertical component of the divergence
			for (int i = 0; i < n; ++i) {
				inner_data [i * m] = (data_z [i * m + 1] - data_z [i * m]) / (pos_z [i * m + 1] - pos_z [i * m]);
				for (int j = 1; j < m - 1; ++j) {
					inner_data [i * m + j] = (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / (pos_z [i * m + j + 1] - pos_z [i * m + j - 1]);
				}
				inner_data [i * m + m - 1] = (data_z [i * m + m - 1] - data_z [i * m + m - 2]) / (pos_z [i * m + m - 1] - pos_z [i * m + m - 2]);
			}
			// Calculate the horizontal component of the divergence
			for (int j = 0; j < m; ++j) {
				inner_data [j] += (data_x [m + j] - data_x [(n - 1) * m + j]) / (pos_x [m + j] - pos_x [(n - 1) * m + j]);
				for (int i = 1; i < n - 1; ++i) {
					inner_data [i * m + j] += (data_x [(i + 1) * m + j] - data_x [(i - 1) * m + j]) / (pos_x [(i + 1) * m + j] - pos_x [(i - 1) * m + j]);
				}
				inner_data [(n - 1) * m + j] += (data_x [j] - data_x [(n - 2) * m + j]) / (pos_x [j] - pos_x [(n - 2) * m + j]);
			}
			return &inner_data [0];
		}
	};
	
	/*!**********************************************************************
	 * \brief Calculates the divergence in transformed space
	 ************************************************************************/
	template <class datatype>
	class transform_div_functor : public functor
	{
	private:
		datatype *pos_x; //!< A datatype pointer to the horizontal position data
		datatype *pos_z; //!< A datatype pointer to the vertical position data
		datatype *data_x; //!< A datatype pointer to the x component of the vector data
		datatype *data_z; //!< A datatype pointer to the z component of the vector data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_pos_x A datatype pointer to the horizontal position data
		 * \param i_pos_z A datatype pointer to the vertical position data
		 * \param i_data_x A datatype pointer to the x component of the vector data
		 * \param i_data_z A datatype pointer to the z component of the vector data
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		transform_div_functor (datatype *i_pos_x, datatype *i_pos_z, datatype *i_data_x, datatype *i_data_z, int i_n, int i_m) : pos_x (i_pos_x), pos_z (i_pos_z), data_x (i_data_x), data_z (i_data_z), n (i_n), m (i_m) {
			inner_data.resize (n * m);
		}

		/*!**********************************************************************
		 * \brief Take the divergence of the data and return a pointer to the first element
		 * 
		 * \return The first element of the 2D divergence array
		 ************************************************************************/
		void *calculate () {
			datatype scalar = acos (-1.0) * 2.0 / (pos_x [m * (n - 1)] - pos_x [0]);
			for (int i = 0; i < n; ++i) {
				inner_data [i * m] = (data_z [i * m + 1] - data_z [i * m]) / (pos_z [i * m + 1] - pos_z [i * m]);
				for (int j = 1; j < m - 1; ++j) {
					inner_data [i * m + j] = (data_z [i * m + j + 1] - data_z [i * m + j - 1]) / (pos_z [i * m + j + 1] - pos_z [i * m + j - 1]);
				}
				inner_data [i * m + m - 1] = (data_z [i * m + m - 1] - data_z [i * m + m - 2]) / (pos_z [i * m + m - 1] - pos_z [i * m + m - 2]);
			}
			for (int j = 0; j < m; ++j) {
				for (int i = 2; i < n; i += 2) {
					inner_data [i * m + j] += -scalar * (i / 2) * data_x [(i + 1) * m + j];
					inner_data [(i + 1) * m + j] += scalar * (i / 2) * data_x [i * m + j];
				}
			}
			return &inner_data [0];
		}
	};
} /* functors */

#endif /* end of include guard: DIV_HPP_A5256DE8 */
