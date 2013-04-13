/*!***********************************************************************
 * \file plan.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PLAN_HPP_S9YPWHOM
#define PLAN_HPP_S9YPWHOM

/*!*******************************************************************
 * \brief The basic functional unit, containing a recipe for execution
 * 
 * An implemented plan class contains the operator and the addresses of 
 * all the relevant data arrays to operate on. Each plan need be 
 * constructed only once and can run any number of times each timestep.
 *********************************************************************/
class plan
{
public:
	virtual ~plan () {}
	/*!*******************************************************************
	 * \brief Operate the plan on the data arrays contained in the class
	 * 
	 * The plan class serves as a wrapper for this function. The user 
	 * specifies the time step and the plan determines how to take care of 
	 * the operation itself.
	 * 
	 * \param timestep A double length of time over which to operate
	 *********************************************************************/
	virtual void execute (double timestep) = 0;
};

#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
