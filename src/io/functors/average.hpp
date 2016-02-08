/*!**********************************************************************
 * \file average.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef AVERAGE_HPP_9C2B188F
#define AVERAGE_HPP_9C2B188F

#include <vector>
#include <memory>

#include "functor.hpp"

namespace functors
{
	/*!**********************************************************************
	 * \brief Averages a two dimensional block of data weighted by an array
	 ************************************************************************/
	template <class datatype>
	class weighted_average_functor : public functor
	{
	private:
		datatype *weight, *data; //!< A datatype pointer to the input data
		std::shared_ptr <functor> func; //!< A shared pointer to another functor object, for nesting
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		datatype inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_weight The datatype pointer to the weighting array
		 * \param i_data The datatype pointer to the data to average
		 ************************************************************************/
		weighted_average_functor (int i_n, int i_m, datatype *i_weight, datatype *i_data) : weight (i_weight), data (i_data), n (i_n), m (i_m) {
		}
		
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_weight The datatype pointer to the weighting array
		 * \param i_func The shared pointer to the functor containing the processed data
		 ************************************************************************/
		weighted_average_functor (int i_n, int i_m, datatype *i_weight, std::shared_ptr <functor> i_func) : weight (i_weight), data ((datatype *) i_func->calculate ()), func (i_func), n (i_n), m (i_m) {
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
			datatype sum = 0.0;
			for (int i = 0; i < m * n; ++i) {
				sum += weight [i];
			}
			inner_data = linalg::dot (m * n, weight, data) / sum;
			// DEBUG ("INNER DATA = " << inner_data);
			// if (inner_data > 1.e10) {
			// 	for (int i = 0; i < m * n; ++i) {
			// 		DEBUG ("WEIGHT IS " << weight[i]);
			// 		DEBUG ("DATUM IS " << data[i]);
			// 	}
			// }
			return &inner_data;
		}
	};

	/*!**********************************************************************
	 * \brief Averages a two dimensional block of data in the horizontal direction
	 ************************************************************************/
	template <class datatype>
	class average_functor : public functor
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		std::shared_ptr <functor> func; //!< A shared pointer to another functor object, for nesting
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A vector of processed data to output
	
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to average
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		average_functor (datatype *i_data, int i_n, int i_m) : data (i_data), n (i_n), m (i_m) {
			inner_data.resize (m);
		}
		
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_func The shared pointer to the functor containing the processed data
		 ************************************************************************/
		average_functor (int i_n, int i_m, std::shared_ptr <functor> i_func) : data ((datatype *) i_func->calculate ()), func (i_func), n (i_n), m (i_m) {
			inner_data.resize (m);
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
			for (int j = 0; j < m; ++j) {
				inner_data [j] = (datatype) 0;
				for (int i = 0; i < n; ++i) {
					inner_data [j] += data [i * m + j];
				}
				inner_data [j] /= (datatype) n;
			}
			return &inner_data [0];
		}
	};
	
	/*!**********************************************************************
	 * \brief Calculate the vertical derivatives of the system for every point
	 ************************************************************************/
	template <class datatype>
	class deriv_functor : public functor
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		const datatype *position_m; //!< A pointer to the vertical position data
		std::vector <datatype> inner_data; //!< A vector of processed data to output
	
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to average
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_position_m The datatype pointer to the vertical position data
		 ************************************************************************/
		deriv_functor (datatype *i_data, int i_n, int i_m, const datatype *i_position_m) : data (i_data), n (i_n), m (i_m), position_m (i_position_m) {
			inner_data.resize (m * n);
		}
	
		/*!**********************************************************************
		 * \brief Calculate the vertical derivatives and return a pointer to the first element
		 * 
		 * \return The first element of the averaged 1D array
		 ************************************************************************/
		void *calculate () {
			for (int j = 1; j < m - 1; ++j) {
				for (int i = 0; i < n; ++i) {
					inner_data [i * m + j] = (data [i * m + j + 1] - data [i * m + j - 1]) / (position_m [j + 1] - position_m [j - 1]);
				}
			}
			return &inner_data [0];
		}
	};

	/*!**********************************************************************
	 * \brief Finds the root-mean-square of data in the horizontal direction
	 ************************************************************************/
	template <class datatype>
	class root_mean_square_functor : public functor
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		std::vector <datatype> inner_data; //!< A datatype vector of processed data to output
	
	public:
		/*!**********************************************************************
		 * \param i_data The datatype pointer to the data to root-mean-square
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 ************************************************************************/
		root_mean_square_functor (datatype *i_data, int i_n, int i_m) : data (i_data), n (i_n), m (i_m) {
			inner_data.resize (m);
		}
	
		/*!**********************************************************************
		 * \brief Take the root-mean-square and return a pointer to the first element
		 * 
		 * \return The first element of the resulting array
		 ************************************************************************/
		void *calculate () {
			for (int j = 0; j < m; ++j) {
				inner_data [j] = (datatype) 0;
				for (int i = 0; i < n; ++i) {
					inner_data [j] += data [i * m + j] * data [i * m + j];
				}
				inner_data [j] = sqrt(inner_data [j] / (datatype) n);
			}
			return &inner_data [0];
		}
	};
	
	/*!**********************************************************************
	 * \brief Calculates the maximum of a two-dimensional array
	 ************************************************************************/
	template <class datatype>
	class max_functor : public functor
	{
	private:
		datatype *data; //!< A datatype pointer to the input data
		std::shared_ptr <functor> func; //!< A shared pointer to another functor object, for nesting
		int n; //!< The integer horizontal extent of the data
		int m; //!< The integer vertical extent of the data
		datatype inner_data; //!< A vector of processed data to output

	public:
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_data The datatype pointer to the data to average
		 ************************************************************************/
		max_functor (int i_n, int i_m, datatype *i_data) : data (i_data), n (i_n), m (i_m) {
		}
		
		/*!**********************************************************************
		 * \param i_n The integer horizontal extent of the data
		 * \param i_m The integer vertical extent of the data
		 * \param i_func The shared pointer to the functor containing the processed data
		 ************************************************************************/
		max_functor (int i_n, int i_m, std::shared_ptr <functor> i_func) : data ((datatype *) i_func->calculate ()), func (i_func), n (i_n), m (i_m) {
		}
		
		/*!**********************************************************************
		 * \brief Find the maximum value of the data and return a pointer to it
		 * 
		 * \return A pointer to the maximum value
		 ************************************************************************/
		void *calculate () {
			if (func) {
				func->calculate ();
			}
			inner_data = 0.0;
			for (int i = 0; i < m * n; ++i) {
				if (data [i] > inner_data) inner_data = data [i];
			}
			return &inner_data;
		}
	};
} /* functors */

#endif /* end of include guard: AVERAGE_HPP_9C2B188F */
