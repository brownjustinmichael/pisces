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
		class implemented_transformer : public plans::transforms::transformer
		{
		private:
			using plans::transforms::transformer::component_flags;
			
			int ldn; //!< The horizontal extent of the data array
			int ldm; //!< The vertical extent of the data array
			grids::variable &data; //!< A pointer to the data to input
			
			std::shared_ptr <plans::plan> forward_horizontal_transform; //!< A shared pointer to the forward horizontal transform
			std::shared_ptr <plans::plan> forward_vertical_transform; //!< A shared pointer to the forward vertical transform
			std::shared_ptr <plans::plan> inverse_horizontal_transform; //!< A shared pointer to the inverse horizontal transform
			std::shared_ptr <plans::plan> inverse_vertical_transform; //!< A shared pointer to the inverse vertical transform
			
		public:
			/*!**********************************************************************
			 * \copydoc transformer::transformer
			 * 
			 * \param i_grid_n The horizontal grid object
			 * \param i_grid_m The vertical grid object
			 * \param i_data A pointer to the data to input and output
			 * \param i_flags Integer flags to describe the setup (e.g. forward_vertical, inverse_horizontal, etc.)
			 * \param i_threads The number of threads to use in the transform
			 ************************************************************************/
			implemented_transformer (grids::variable &i_data, int *i_element_flags, int *i_component_flags, int i_flags = forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal, int i_threads = 1) : 
			plans::transforms::transformer (i_element_flags, i_component_flags), 
			ldn (i_data.get_grid (0).get_ld ()), 
			ldm (i_data.get_grid (1).get_ld ()), 
			data (i_data) {
				// For each direction, check the flags to see which transforms to add and do so
				TRACE("Initializing...");
				if (i_flags & forward_vertical) {
					forward_vertical_transform = std::shared_ptr <plans::plan> (new plans::transforms::vertical (data, real_spectral, spectral_spectral, 0x00, i_threads));
				}
				if (i_flags & inverse_vertical) {
					inverse_vertical_transform = std::shared_ptr <plans::plan> (new plans::transforms::vertical (data, spectral_spectral, real_spectral, inverse, i_threads));
				}
				if (i_flags & forward_horizontal) {
					forward_horizontal_transform = std::shared_ptr <plans::plan> (new plans::transforms::horizontal (data, real_real, real_spectral, 0x00, i_threads));
				}
				if (i_flags & inverse_horizontal) {
					inverse_horizontal_transform = std::shared_ptr <plans::plan> (new plans::transforms::horizontal (data, real_spectral, real_real, inverse, i_threads));
				}
			}
			
			virtual ~implemented_transformer () {}
			
			/*!**********************************************************************
			 * \copydoc transformer::update
			 ************************************************************************/
			void update () {
				if ((data.component_flags & grids::updated)) {
					DEBUG ("Skipping");
					return;
				}
				int state = data.last_update;
				DEBUG ("Updating from " << state);
				if (state == real_real) {
					forward_horizontal_transform->execute ();
					forward_vertical_transform->execute ();
				} else if (state == real_spectral) {
					inverse_horizontal_transform->execute ();
					forward_vertical_transform->execute ();
				} else if (state == spectral_spectral) {
					inverse_vertical_transform->execute ();
					inverse_horizontal_transform->execute ();
				}
				data.component_flags |= grids::updated;
				data.state++;
			}
		};
	} /* transforms */
} /* plans */

#endif /* end of include guard: IMPLEMENTED_TRANSFORMER_HPP_8A68B0E1 */
