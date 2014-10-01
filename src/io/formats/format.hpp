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
	
	class data_grid
	{
	private:
		int n_dims;
		
	public:
		std::vector <size_t> ns;
		std::vector <size_t> maxes;
		std::vector <size_t> offsets;
		
		data_grid (int initial_record = 0, int i_n_dims = 0, size_t *i_ns = NULL, size_t *i_maxes = NULL, size_t *i_offsets = NULL) {
			add_dimension (1, 1, initial_record);
			for (int i = 0; i < i_n_dims; ++i) {
				add_dimension (i_ns [i], i_maxes [i], i_offsets [i]);
			}
		}
		
		virtual ~data_grid () {}
		
		int get_n_dims () const {
			return n_dims - 1;
		}
		
		size_t get_n (int i) const {
			return ns [i + 1];
		}
		
		size_t get_max (int i) const {
			return maxes [i + 1];
		}
		
		size_t get_offset (int i) const {
			return offsets [i + 1];
		}
		
		static data_grid one_d (int n, int max = 0, int offset = 0) {
			data_grid out;
			out.add_dimension ((size_t) n, (size_t) max, (size_t) offset);
			return out;
		}
		
		static data_grid two_d (int n, int m, int n_max = 0, int m_max = 0, int n_offset = 0, int m_offset = 0) {
			data_grid out;
			out.add_dimension ((size_t) n, (size_t) n_max, (size_t) n_offset);
			out.add_dimension ((size_t) m, (size_t) m_max, (size_t) m_offset);
			return out;
		}
		
		static data_grid three_d (int n, int m, int l, int n_max = 0, int m_max = 0, int l_max = 0, int n_offset = 0, int m_offset = 0, int l_offset = 0) {
			data_grid out;
			out.add_dimension ((size_t) n, (size_t) n_max, (size_t) n_offset);
			out.add_dimension ((size_t) m, (size_t) m_max, (size_t) m_offset);
			out.add_dimension ((size_t) l, (size_t) l_max, (size_t) l_offset);
			return out;
		}
	private:
		void add_dimension (size_t n, size_t max, size_t offset) {
			n_dims += 1;
			ns.push_back (n);
			maxes.push_back (max > 0 ? max : n);
			offsets.push_back (offset);
		}
	};
} /* io */

#endif /* end of include guard: FORMAT_HPP_6ADF5F85 */
