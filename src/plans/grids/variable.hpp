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
		int states;
		// TODO This won't work. We need a way to keep references or entire variables
		std::vector <variable <datatype> *> inner;
		std::vector <int> ops;
		std::vector <int> inner_states;
		int ld;
		int total;

	public:
		int state = 0;
		int last_update = 0;
		int component_flags;
		int &element_flags;
		std::string name;

		enum operators {
			add = 0x01,
			sub = 0x02,
			mul = 0x03,
			div = 0x04
		};

		variable (grids::grid <datatype> &i_grid_m, int &i_element_flags, std::string i_name = "", int i_states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags), name (i_name) {
			DEBUG ("Initializing " << name << " with " << i_states);
			grids.push_back (&i_grid_m);
			total = i_grid_m.get_n ();
			DEBUG ("TOTAL " << total);
			data.resize (total * dimensions * i_states, 0.0);
			ld = 0;
			states = i_states;
		}

		variable (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, int &i_element_flags, std::string i_name = "", int i_states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags), name (i_name) {
			DEBUG ("Initializing " << name << " with " << i_states);
			grids.push_back (&i_grid_n);
			grids.push_back (&i_grid_m);
			total = i_grid_m.get_n () * i_grid_n.get_ld ();
			DEBUG ("TOTAL " << total);
			data.resize (total * i_states * dimensions, 0.0);
			ld = i_grid_m.get_n ();
			states = i_states;
		}

		variable (int n, grids::grid <datatype> **i_grids, int &i_element_flags, std::string i_name = "", int i_states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags), name (i_name) {
			DEBUG ("Initializing " << name << " with " << i_states);
			total = n >= 1 ? 1 : 0;
			for (int i = 0; i < n; ++i)
			{
				grids.push_back (i_grids [i]);
				total *= i_grids [i]->get_ld ();
			}
			DEBUG ("TOTAL " << total);
			data.resize (total * i_states * dimensions, 0.0);
			ld = i_grids [n - 1]->get_ld ();
			states = i_states;
		}

		virtual ~variable () {}

		datatype &operator[] (int index) {
			component_flags &= ~updated;
			last_update = 0;
			state++;
			return data [index];
		}

		datatype *ptr (int i_state = 0) {
			if (i_state >= states) {
				ERROR ("State " << i_state << " not initialized.");
				throw 501;
			}
			return &data [i_state * total * dimensions];
		}

		std::string get_name ();

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

		const int get_states () {
			return states;
		}


		int get_ld () {
			return ld;
		}

		void add_var (variable <datatype> &data, int op) {
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
			this->update ();
			return *this;
		}

		static void store_var (std::shared_ptr <variable <datatype>> &other) {
			tmps.push_back (other);
		}

		static void update_tmps () {
			for (auto iter = tmps.begin (); iter != tmps.end (); ++iter)
			{
				(*iter)->update ();
			}
		}

		variable <datatype> &operator+ (datatype other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			std::shared_ptr <variable <datatype>> uni_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, std::to_string (other), 1)));
			linalg::copy (this->size (), &other, uni_var->ptr (), 0);
			new_var->add_var (*this, add);
			new_var->add_var (*uni_var, add);
			store_var (new_var);
			store_var (uni_var);
			new_var->update ();

			return *new_var;
		}

		variable <datatype> &operator+ (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			new_var->add_var (*this, add);
			new_var->add_var (other, add);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		variable <datatype> &operator- (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			new_var->add_var (*this, add);
			new_var->add_var (other, sub);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		variable <datatype> &operator* (datatype other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			std::shared_ptr <variable <datatype>> uni_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, std::to_string (other), 1)));
			linalg::copy (this->size (), &other, uni_var->ptr (), 0);
			new_var->add_var (*this, add);
			new_var->add_var (*uni_var, mul);
			store_var (new_var);
			store_var (uni_var);
			new_var->update ();

			return *new_var;
		}

		variable <datatype> &operator* (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var;
			if (this->shape () > other.shape ()) {
				new_var.reset (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, "", 1));
			} else {
				new_var.reset (new variable <datatype> (other.shape (), other.get_grids (), this->element_flags, "", 1));
			}			new_var->add_var (*this, add);
			new_var->add_var (other, mul);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		variable <datatype> &operator/ (variable <datatype> &other) {
			DEBUG ("THE PREVIOUS SHAPE IS " << this->shape ());
			std::shared_ptr <variable <datatype>> new_var;
			if (this->shape () > other.shape ()) {
				new_var.reset (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags, "", 1));
			} else {
				new_var.reset (new variable <datatype> (other.shape (), other.get_grids (), this->element_flags, "", 1));
			}
			new_var->add_var (*this, add);
			new_var->add_var (other, div);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}
	};

	template <class datatype>
	variable <datatype> &operator+ (datatype first, variable <datatype> &other) {
		return other + first;
	}

	template <class datatype>
	variable <datatype> &operator- (datatype first, variable <datatype> &other) {
		return other * (-1.) + first;
	}

	template <class datatype>
	variable <datatype> &operator/ (datatype first, variable <datatype> &other) {
		std::shared_ptr <variable <datatype>> uni_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (other.shape (), other.get_grids (), other.element_flags, std::to_string (first))));
		linalg::copy (other.size (), &first, uni_var->ptr (), 0);
		variable <datatype>::store_var (uni_var);

		return *uni_var / other;
	}
}

#endif /* end of include guard: VARIABLE_H__ */
