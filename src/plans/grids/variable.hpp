/*!***********************************************************************
 * \file variable.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VARIABLE_H__
#define VARIABLE_H__

#include "grid.hpp"

namespace grids
{
	enum variable_flags
	{
		updated = 0x80000000
	};

	template <class datatype>
	class variable
	{
	protected:
		static std::vector <std::shared_ptr <variable <datatype>>> tmps;

		std::vector <grids::grid <datatype> *> grids;
		int dimensions;
		std::vector <datatype> data;
		std::vector <datatype *> vars;
		// TODO This won't work. We need a way to keep references or entire variables
		std::vector <variable <datatype> *> inner;
		std::vector <int> ops;
		std::vector <int> inner_states;
		int ld;
		int total;

	public:
		int state = 0;
		int last_update;
		int component_flags;
		int &element_flags;

		enum operators {
			add = 0x01,
			sub = 0x02,
			mul = 0x03,
			div = 0x04
		};

		variable (grids::grid <datatype> &i_grid_m, int &i_element_flags, int states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags) {
			grids.push_back (&i_grid_m);
			total = i_grid_m.get_n ();
			data.resize (total * dimensions * states, 0.0);
			ld = 0;
		}

		variable (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, int &i_element_flags, int states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags) {
			grids.push_back (&i_grid_n);
			grids.push_back (&i_grid_m);
			data.resize (i_grid_m.get_n () * i_grid_n.get_ld () * states * dimensions, 0.0);
			ld = i_grid_m.get_n ();
			total = i_grid_m.get_n () * i_grid_n.get_ld ();
		}

		variable (int n, grids::grid <datatype> **i_grids, int &i_element_flags, int states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags) {
			total = n > 1 ? 1 : 0;
			for (int i = 0; i < n; ++i)
			{
				grids.push_back (i_grids [i]);
				total *= i_grids [i]->get_ld ();
			}
			data.resize (total * states * dimensions, 0.0);
			ld = i_grids [n - 1]->get_ld ();
		}

		virtual ~variable () {}

		datatype *ptr (int state = 0) {
			return &data [state * total * dimensions];
		}

		const int size () const {
			return total * dimensions;
		}

		const int shape () const {
			return grids.size ();
		}

		grid <datatype> &get_grid (int n) {
			return *grids [n];
		}

		grid <datatype>** get_grids () {
			return &grids [0];
		}

		const int dims () {
			return dimensions;
		}

		int get_ld () {
			return ld;
		}

		void add_var (variable <datatype> &data, int op) {
			DEBUG ("ADDING " << data.ptr ());
			inner.push_back (&data);
			vars.push_back (data.ptr ());
			ops.push_back (op);
			inner_states.push_back (-1);
		}

		void reset_vars () {
			inner.resize (0);
			vars.resize (0);
			ops.resize (0);
			inner_states.resize (0);
		}

		int update ();

		variable <datatype> &operator== (variable <datatype> &other) {
			this->reset_vars ();
			this->add_var (other, add);
			return *this;
		}

		variable <datatype> &operator* (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags)));
			new_var->add_var (*this, add);
			new_var->add_var (other, mul);
			tmps.push_back (new_var);

			return *new_var;
		}

		variable <datatype> &operator/ (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags)));
			new_var->add_var (*this, add);
			new_var->add_var (other, div);
			tmps.push_back (new_var);

			return *new_var;
		}
	};
}

#endif /* end of include guard: VARIABLE_H__ */
