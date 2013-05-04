/*!***********************************************************************
 * \file timestep.hpp
 * src
 * 
 * Created by Justin Brown on 2013-04-26.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plan.hpp"
#include "../config.hpp"

namespace bases
{
	class calculate_timestep : public plan
	{
	public:
		calculate_timestep (double& i_timestep) : timestep (i_timestep) {
			MTRACE ("Instantiated.");
		}
		virtual ~calculate_timestep () {}
		
		virtual void execute () = 0;

	protected:
		double& timestep;
	};
} /* bases */

class constant_timestep : public bases::calculate_timestep
{
public:
	constant_timestep (double i_initial_timestep, double& i_timestep) : calculate_timestep (i_timestep) {
		initial_timestep = i_initial_timestep;
	}
	
	virtual ~constant_timestep () {}
	
	virtual void execute () {
		timestep = initial_timestep;
	}

private:
	double initial_timestep;
};