/*!***********************************************************************
 * \file variable.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "variable.hpp"

namespace grids
{
	std::vector <std::shared_ptr <variable>> grids::variable::tmps;

	int variable::update () {
		bool needs_update = false;
		int new_state;
		for (int i = 0; i < (int) vars.size (); ++i)
		{
			new_state = inner [i]->update ();
			needs_update = needs_update || (new_state != inner_states [i]);
			inner_states [i] = new_state;
		}

		if (!needs_update) {
			DEBUG ("No update needed");
			return state;
		}

		DEBUG ("Update attempt on " << this->get_name ());
		linalg::scale (this->size (), 0.0, this->ptr ());
		double *tmp, *data = this->ptr ();
		for (int i = 0; i < (int) vars.size (); ++i)
		{
			DEBUG ("Incorporating " << inner [i]->get_name () << " at " << inner [i]->ptr ());
			tmp = inner [i]->ptr ();
			int n = total / ld;
			int tld = inner [i]->get_ld ();
			if (ops [i] == add) {
				for (int i = 0; i < n; ++i)
				{
					linalg::add_scaled (ld, tmp + i * tld, data + i * ld);
				}
			} else if (ops [i] == sub) {
				for (int i = 0; i < n; ++i)
				{
					linalg::add_scaled (ld, -1.0, tmp + i * tld, data + i * ld);
				}
			} else if (ops [i] == mul) {
				#pragma omp parallel for
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < ld; ++j)
					{
						data [i * ld + j] *= tmp [i * tld + j];
					}
				}
			} else if (ops [i] == div) {
				#pragma omp parallel for
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < ld; ++j)
					{
						data [i * ld + j] /= tmp [i * tld + j];
					}
				}
			}
		}
		component_flags &= ~updated;
		last_update = 0;
		state++;
		return state;
	}

	std::string variable::get_name () {
		std::string return_val = name;
		bool inner_bool = false;
		for (int i = 0; i < (int) this->inner.size (); ++i)
		{
		 	inner_bool = true;
		 	if (i != 0) {
		 		switch (this->ops [i]) {
		 			case add:
		 				return_val += "+";
		 				break;
	 				case sub:
	 					return_val += "-";
	 					break;
	 				case mul:
	 					return_val += "*";
	 					break;
	 				case div:
	 					return_val += "/";
	 					break;
		 		}
		 	} else {
		 		return_val += "(";
		 	}
		 	return_val += this->inner [i]->get_name ();
		}
		if (inner_bool) return_val += ")";
		return return_val;
	}

	variable &operator+ (double first, variable &other) {
		return other + first;
	}

	variable &operator- (double first, variable &other) {
		return other * (-1.) + first;
	}

	variable &operator/ (double first, variable &other) {
		std::shared_ptr <variable> uni_var (std::shared_ptr <variable> (new variable (other.shape (), other.get_grids (), other.element_flags, std::to_string ((long double) first))));
		linalg::copy (other.size (), &first, uni_var->ptr (), 0);
		variable::store_var (uni_var);

		return *uni_var / other;
	}
}