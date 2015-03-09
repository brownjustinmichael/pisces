/*!**********************************************************************
 * \file format.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FORMAT_HPP_6ADF5F85
#define FORMAT_HPP_6ADF5F85

#include <vector>

#include "versions/version.hpp"

namespace io
{
	/*!**********************************************************************
	 * \brief A set of io flags to be used with the input/output classes specifying the file type
	 ************************************************************************/
	enum io_flags {
		read_file = 0,
		replace_file = 1,
		append_file = 2
	};
	
	/*!**********************************************************************
	 * \brief A set of dim flags to be used with the input/output classes specifying the dimensionality of the data
	 ************************************************************************/
	enum dim_flags {
		all_d = -0x1,
		scalar = 0x0,
		one_d = 0x1,
		m_profile = 0x2,
		two_d = 0x3,
		l_profile = 0x4,
		three_d = 0x7
	};
	
	/*!**********************************************************************
	 * \brief An object that contains the relevant sizes of dimensions for the data
	 * 
	 * This object is designed for use with input/output classes in order to keep track of the dimensionality of the data. The offsets here are defined to allow for each element of a spectral element code to only write to a small part of an output file, allowing for output files to be stacked in a sensible way.
	 ************************************************************************/
	class data_grid
	{
	private:
		int n_dims; //!< The number of dimensions in the data
		
	public:
		std::vector <size_t> ns; //!< The actual extent of the data in each dimension
		std::vector <size_t> maxes; //!< The full extent of the arrays in the file for each dimension
		std::vector <size_t> offsets; //!< The offset of the actual data and the file arrays for each dimension
		
		/*!**********************************************************************
		 * \param initial_record The initial temporal record number of the data grid 
		 * \param i_n_dims The number of dimensions to output
		 * \param i_ns An array of the extents of each dimension to output
		 * \param i_maxes The maximum extent of each dimension in the file
		 * \param i_offsets The offsets of the beginning point of each dimension in the file
		 ************************************************************************/
		data_grid (int initial_record = 0, int i_n_dims = 0, size_t *i_ns = NULL, size_t *i_maxes = NULL, size_t *i_offsets = NULL) {
			n_dims = 0;
			add_dimension (1, 1, initial_record);
			for (int i = 0; i < i_n_dims; ++i) {
				add_dimension (i_ns [i], i_maxes [i], i_offsets [i]);
			}
		}
		
		virtual ~data_grid () {}
		
		/*!**********************************************************************
		 * \brief Get the version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.0.1.0");
			return version;
		}
		
		/*!**********************************************************************
		 * \return The number of dimensions
		 ************************************************************************/
		int get_n_dims () const {
			return n_dims - 1;
		}
		
		/*!**********************************************************************
		 * \param i The dimension to index
		 * 
		 * \return The extent of the data in dimension i
		 ************************************************************************/
		size_t get_n (int i) const {
			return ns [i + 1];
		}
		
		/*!**********************************************************************
		 * \param i The dimension to index
		 * 
		 * \return The maximum extent of the file array in dimension i
		 ************************************************************************/
		size_t get_max (int i) const {
			return maxes [i + 1];
		}
		
		/*!**********************************************************************
		 * \param i The dimension to index
		 * 
		 * \return The offset of the data in the file
		 ************************************************************************/
		size_t get_offset (int i) const {
			return offsets [i + 1];
		}
		
		/*!**********************************************************************
		 * \brief A static constructor to create a 1D grid
		 * 
		 * \param n The extent of the data
		 * \param max The max extent of the array
		 * \param offset The offset of the data in the array
		 * 
		 * \return The constructed 1D grid
		 ************************************************************************/
		static data_grid one_d (int n, int max = 0, int offset = 0) {
			data_grid out;
			out.add_dimension ((size_t) n, (size_t) max, (size_t) offset);
			return out;
		}
		
		/*!**********************************************************************
		 * \brief A static constructor to create a 2D grid
		 * 
		 * \param n The extent of the data in dimension 1
		 * \param m The extent of the data in dimension 2
		 * \param n_max The max extent of the array in dimension 1
		 * \param m_max The max extent of the array in dimension 2
		 * \param n_offset The offset of the data in the array in dimension 1
		 * \param m_offset The offset of the data in the array in dimension 2
		 * 
		 * \return The constructed 2D grid
		 ************************************************************************/
		static data_grid two_d (int n, int m, int n_max = 0, int m_max = 0, int n_offset = 0, int m_offset = 0) {
			data_grid out;
			out.add_dimension ((size_t) n, (size_t) n_max, (size_t) n_offset);
			out.add_dimension ((size_t) m, (size_t) m_max, (size_t) m_offset);
			return out;
		}
		
		/*!**********************************************************************
		 * \brief A static constructor to create a 3D grid
		 * 
		 * \param n The extent of the data in dimension 1
		 * \param m The extent of the data in dimension 2
		 * \param l The extent of the data in dimension 3
		 * \param n_max The max extent of the array in dimension 1
		 * \param m_max The max extent of the array in dimension 2
		 * \param l_max The max extent of the array in dimension 3
		 * \param n_offset The offset of the data in the array in dimension 1
		 * \param m_offset The offset of the data in the array in dimension 2
		 * \param l_offset The offset of the data in the array in dimension 3
		 * 
		 * \return The constructed 3D grid
		 ************************************************************************/
		static data_grid three_d (int n, int m, int l, int n_max = 0, int m_max = 0, int l_max = 0, int n_offset = 0, int m_offset = 0, int l_offset = 0) {
			data_grid out;
			out.add_dimension ((size_t) n, (size_t) n_max, (size_t) n_offset);
			out.add_dimension ((size_t) m, (size_t) m_max, (size_t) m_offset);
			out.add_dimension ((size_t) l, (size_t) l_max, (size_t) l_offset);
			return out;
		}
		
	private:
		/*!**********************************************************************
		 * \brief Add a dimension to the data_grid
		 * 
		 * \param n The extent of the new dimension
		 * \param max The max array length of the new dimension
		 * \param offset The offset of the data in the array
		 * 
		 * This is declared private to prevent changing the array structure after definition.
		 ************************************************************************/
		void add_dimension (size_t n, size_t max, size_t offset) {
			n_dims += 1;
			ns.push_back (n);
			maxes.push_back (max > 0 ? max : n);
			offsets.push_back (offset);
		}
	};
} /* io */

#endif /* end of include guard: FORMAT_HPP_6ADF5F85 */
