/*!**********************************************************************
 * \file implemented_transformer.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-07.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IMPLEMENTED_TRANSFORMER_HPP_8A68B0E1
#define IMPLEMENTED_TRANSFORMER_HPP_8A68B0E1

#include "plans/plan.hpp"

#include "transformer.hpp"
#include "transform.hpp"

namespace plans
{
	template <class datatype>
	class implemented_transformer : public plans::transformer <datatype>
	{
	public:
		implemented_transformer (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags, int i_threads) : 
		plans::transformer <datatype> (element_flags, component_flags),
		ldn (i_grid_n.get_ld ()),
		ldm (i_grid_m.get_ld ()),
		data_in (i_data_in),
		data_out (i_data_out ? i_data_out : i_data_in) {
			data.resize (ldn * ldm, 0.0);
			if (i_flags & forward_vertical) {
				forward_vertical_transform = std::shared_ptr <plans::plan <datatype>> (new plans::vertical_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, element_flags, &internal_state, i_threads));
			}
			if (i_flags & inverse_vertical) {
				if (forward_vertical_transform) {
					inverse_vertical_transform = forward_vertical_transform;
				} else {
					inverse_vertical_transform = std::shared_ptr <plans::plan <datatype>> (new plans::vertical_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, element_flags, &internal_state, i_threads));
				}
			}
			if (i_flags & forward_horizontal) {
				forward_horizontal_transform = std::shared_ptr <plans::plan <datatype>> (new plans::horizontal_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, element_flags, &internal_state, i_threads));
			}
			if (i_flags & inverse_horizontal) {
				inverse_horizontal_transform = std::shared_ptr <plans::plan <datatype>> (new plans::horizontal_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, element_flags, &internal_state, i_threads));
			}
		}
		
		virtual ~implemented_transformer () {}
	
		void _transform (int flags) {
			if (flags & forward_horizontal) {
				if (!(internal_state & transformed_horizontal) && forward_horizontal_transform) {
					forward_horizontal_transform->execute ();
				}
			}
			if (flags & forward_vertical) {
				if (!(internal_state & transformed_vertical) && forward_vertical_transform) {
					forward_vertical_transform->execute ();
				}
			}
			if (flags & inverse_horizontal) {
				if ((internal_state & transformed_horizontal) && inverse_horizontal_transform) {
					inverse_horizontal_transform->execute ();
				}
			}
			if (flags & inverse_vertical) {
				if ((internal_state & transformed_vertical) && inverse_vertical_transform) {
					inverse_vertical_transform->execute ();
				}
			}
		}
	
		void write () {
			linalg::matrix_copy (ldm, ldn, data_in, &data [0]);
		}

		void read () {
			linalg::matrix_copy (ldm, ldn, &data [0], data_out);
		}
		
		virtual datatype *get_data () {
			return &data [0];
		}
	
	private:
		int ldn, ldm;
		datatype *data_in, *data_out;
		std::vector <datatype> data;
		
		std::shared_ptr <plans::plan <datatype>> forward_horizontal_transform;
		std::shared_ptr <plans::plan <datatype>> forward_vertical_transform;
		std::shared_ptr <plans::plan <datatype>> inverse_horizontal_transform;
		std::shared_ptr <plans::plan <datatype>> inverse_vertical_transform;
		
		using plans::transformer <datatype>::internal_state;
		using plans::transformer <datatype>::component_flags;
	};
} /* plans */

#endif /* end of include guard: IMPLEMENTED_TRANSFORMER_HPP_8A68B0E1 */
