/*!***********************************************************************
 * \file timestep_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-05-03.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TIMESTEP_ONE_D_HPP_XRLQJFKX
#define TIMESTEP_ONE_D_HPP_XRLQJFKX

#include "../bases/timestep.hpp"

namespace one_d
{
	class advection_timestep : public bases::calculate_timestep
	{
	public:
		calculate_timestep (double i_cfl, double i_advection_coefficient, double& i_timestep, double* i_position, double* i_velocity);
		virtual ~calculate_timestep ();
	
	private:
		/* data */
	};
} /* one_d */

#endif /* end of include guard: TIMESTEP_ONE_D_HPP_XRLQJFKX */
