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
	namespace transforms
	{
		/*!**********************************************************************
		 * \brief An implemented transformer object in 2D
		 ************************************************************************/
		template <class datatype>
		class implemented_transformer : public plans::transforms::transformer <datatype>
		{
		private:
			using plans::transforms::transformer <datatype>::internal_state;
			using plans::transforms::transformer <datatype>::component_flags;
			
			int ldn; //!< The horizontal extent of the data array
			int ldm; //!< The vertical extent of the data array
			datatype *data_in; //!< A pointer to the data to input
			datatype *data_out; //!< A pointer to the data to output
			std::vector <datatype> data; //!< A vector containing the inner data
			
			std::shared_ptr <plans::plan <datatype>> forward_horizontal_transform; //!< A shared pointer to the forward horizontal transform
			std::shared_ptr <plans::plan <datatype>> forward_vertical_transform; //!< A shared pointer to the forward vertical transform
			std::shared_ptr <plans::plan <datatype>> inverse_horizontal_transform; //!< A shared pointer to the inverse horizontal transform
			std::shared_ptr <plans::plan <datatype>> inverse_vertical_transform; //!< A shared pointer to the inverse vertical transform
			
		public:
			/*!**********************************************************************
			 * \copydoc transformer::transformer
			 * 
			 * \param i_grid_n The horizontal grid object
			 * \param i_grid_m The vertical grid object
			 * \param i_data_in A pointer to the data to input
			 * \param i_data_out A pointer to the data to output
			 * \param i_flags Integer flags to describe the setup (e.g. forward_vertical, inverse_horizontal, etc.)
			 * \param i_threads The number of threads to use in the transform
			 ************************************************************************/
			implemented_transformer (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads) : plans::transforms::transformer <datatype> (i_element_flags, i_component_flags), ldn (i_grid_n.get_ld ()), ldm (i_grid_m.get_ld ()), data_in (i_data_in), data_out (i_data_out ? i_data_out : i_data_in) {
				data.resize (ldn * ldm, 0.0);
				// For each direction, check the flags to see which transforms to add and do so
				if (i_flags & forward_vertical) {
					forward_vertical_transform = std::shared_ptr <plans::plan <datatype>> (new plans::transforms::vertical <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, i_element_flags, &internal_state, i_threads));
				}
				if (i_flags & inverse_vertical) {
					if (forward_vertical_transform) {
						inverse_vertical_transform = forward_vertical_transform;
					} else {
						inverse_vertical_transform = std::shared_ptr <plans::plan <datatype>> (new plans::transforms::vertical <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, i_element_flags, &internal_state, i_threads));
					}
				}
				if (i_flags & forward_horizontal) {
					forward_horizontal_transform = std::shared_ptr <plans::plan <datatype>> (new plans::transforms::horizontal <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, i_element_flags, &internal_state, i_threads));
				}
				if (i_flags & inverse_horizontal) {
					inverse_horizontal_transform = std::shared_ptr <plans::plan <datatype>> (new plans::transforms::horizontal <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, i_element_flags, &internal_state, i_threads));
				}
			}
			
			virtual ~implemented_transformer () {}
			
			/*!**********************************************************************
			 * \copydoc transformer::get_data
			 ************************************************************************/
			virtual datatype *get_data () {
				return &data [0];
			}
			
		protected:
			/*!**********************************************************************
			 * \copydoc transformer::_transform
			 ************************************************************************/
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
			
			/*!**********************************************************************
			 * \copydoc transformer::write
			 ************************************************************************/
			void write () {
				linalg::matrix_copy (ldm, ldn, data_in, &data [0]);
			}
			
			/*!**********************************************************************
			 * \copydoc transformer::read
			 ************************************************************************/
			void read () {
				linalg::matrix_copy (ldm, ldn, &data [0], data_out);
			}
		};
	} /* transforms */
} /* plans */

#endif /* end of include guard: IMPLEMENTED_TRANSFORMER_HPP_8A68B0E1 */
