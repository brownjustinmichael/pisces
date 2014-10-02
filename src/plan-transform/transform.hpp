/*!***********************************************************************
 * \file bases/transform.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2014-02-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_HPP_S8KRHTX3
#define TRANSFORM_HPP_S8KRHTX3
	
#include "plans/grid.hpp"

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

namespace bases
{
	/*!*******************************************************************
	 * \brief A class designed to track and implement the transformations of a particular dataset
	 * 
	 * Note that the design of the element class expects that calling only the transform does not change the dataset. The transformed dataset must first be read back into the original for the transform to take effect.
	 *********************************************************************/
	template <class datatype>
	class master_transform
	{
	protected:
		int *element_flags; //!< A pointer to the flags describing the global state of the element
		int *component_flags; //!< A pointer to the flags describing the state of the local variable

	public:
		/*!**********************************************************************
		 * \param i_element_flags A pointer to the flags describing the global state of the element
		 * \param i_component_flags A pointer to the flags describing the state of the local variable
		 ************************************************************************/
		master_transform (int *i_element_flags, int *i_component_flags) :
		element_flags (i_component_flags),
		component_flags (i_component_flags) {}
		
		virtual ~master_transform () {}
		
		/*!*******************************************************************
		 * \brief Transform the dataset according to the given flags
		 * 
		 * \param flags Binary flags describing the transform, from transform_flags
		 * 
		 * By default, this method will write the data into the transform class, perform the transform, and read the data back out into the element. This is chosen in order to allow for GPU usage in the future.
		 *********************************************************************/
		virtual void transform (int flags = 0x00) {
			if (flags & read_before) {
				read ();
			}
			if (!(flags & no_write)) {
				write ();
			}
			_transform (flags);
			if (!(flags & no_read)) {
				read ();
			}
		}
		
		virtual datatype *get_data () {
			return NULL;
		}
		
	protected:
		/*!**********************************************************************
		 * \brief Write the dataset to the tranform class
		 * 
		 * This must be overwritten in a subclass.
		 ************************************************************************/
		virtual void write () = 0;
		
		/*!**********************************************************************
		 * \brief Read the dataset from the transform class
		 * 
		 * This must be overwritten in a subclass.
		 ************************************************************************/
		virtual void read () = 0;
		
		/*!**********************************************************************
		 * \brief Transform the dataset
		 * 
		 * \params flags Binary flags that determine the type of transform (forward_vertical, forward_horizontal, inverse_vertical, inverse_horizontal)
		 * 
		 * This method contains the implementation of the transform, which must be overwritten in the subclasses.
		 ************************************************************************/
		virtual void _transform (int flags) = 0;
	};
} /* bases */

#endif /* end of include guard: TRANSFORM_HPP_S8KRHTX3 */
