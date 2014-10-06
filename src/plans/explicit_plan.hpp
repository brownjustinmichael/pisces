/*!**********************************************************************
 * \file explicit_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXPLICIT_PLAN_HPP_5171CAFD
#define EXPLICIT_PLAN_HPP_5171CAFD

#include "grid.hpp"
#include "plan.hpp"

namespace plans
{
   /*!*******************************************************************
    * \brief A subclass of plan, specific to explicit methods
    * 
    * These plans take input and produce output.
    *********************************************************************/
   template <class datatype>
   class explicit_plan : public plan <datatype>
   {
   public:
   	explicit_plan (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) :
   	plans::plan <datatype> (i_element_flags, i_component_flags),
   	n (i_grid_n.get_n ()),
   	ldn (i_grid_n.get_ld ()),
   	m (i_grid_m.get_n ()),
   	grid_n (i_grid_n),
   	grid_m (i_grid_m),
   	data_in (i_data_in),
   	data_out (i_data_out ? i_data_out : i_data_in) {}

   	virtual ~explicit_plan () {}

   	/*!*******************************************************************
   	 * \copydoc plans::plan::execute ()
   	 *********************************************************************/
   	virtual void execute () = 0;
	
	class factory
	{
	public:
		virtual ~factory () {}
	
		virtual std::shared_ptr <plans::plan <datatype>> instance (plans::grid <datatype> **grids, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const = 0;
	};

   	int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
   	int ldn;
   	int m;
   	plans::grid <datatype> &grid_n, &grid_m;
   	datatype* data_in; //!< A datatype pointer to the input data
   	datatype* data_out; //!< A datatype pointer to the output data
   };
} /* plans */

#endif /* end of include guard: EXPLICIT_PLAN_HPP_5171CAFD */
