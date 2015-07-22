/*!***********************************************************************
 * \file transformer.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2014-02-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_HPP_S8KRHTX3
#define TRANSFORM_HPP_S8KRHTX3

#include "versions/version.hpp"
#include "plans/grids/grid.hpp"

namespace plans
{
	namespace transforms
	{
		/*!**********************************************************************
		 * \brief A set of flags to be used with setting up transforms
		 ************************************************************************/
		enum transform_flags {
			do_not_transform = 0x00,
			forward_horizontal = 0x01,
			forward_vertical = 0x02,
			inverse_horizontal = 0x04,
			inverse_vertical = 0x08,
			ignore_m = 0x10,
			inverse = 0x20,
			no_write = 0x40,
			no_read = 0x80,
			read_before = 0x100
		};
		
		/*!*******************************************************************
		 * \brief A class designed to track and implement the transformations of a particular dataset
		 * 
		 * Note that the design of the element class expects that calling only the transform does not change the dataset. The transformed dataset must first be read back into the original for the transform to take effect.
		 *********************************************************************/
		template <class datatype>
		class transformer
		{
		protected:
			int *element_flags; //!< A pointer to the flags describing the global state of the element
			int *component_flags; //!< A pointer to the flags describing the state of the local variable

		public:
			/*!**********************************************************************
			 * \param i_element_flags A pointer to the flags describing the global state of the element
			 * \param i_component_flags A pointer to the flags describing the state of the local variable
			 ************************************************************************/
			transformer (int *i_element_flags, int *i_component_flags) : element_flags (i_component_flags), component_flags (i_component_flags) {
			}
		
			virtual ~transformer () {}
		
			/*!**********************************************************************
			 * \return The version of the class
			 ************************************************************************/
			static versions::version& version () {
				static versions::version version ("1.0.1.0");
				return version;
			}

			virtual void update () = 0;
		};
	} /* transforms */
} /* plans */

#endif /* end of include guard: TRANSFORM_HPP_S8KRHTX3 */
