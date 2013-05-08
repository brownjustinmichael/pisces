/*!***********************************************************************
 * \file timestep.hpp
 * src
 * 
 * Created by Justin Brown on 2013-04-26.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "plan.hpp"
#include "../config.hpp"
#include "solver.hpp"

namespace bases
{
	class calculate_timestep : public plan
	{
	public:
		calculate_timestep (double& i_timestep) : timestep (i_timestep) {
			previous_timestep = 0.0;
			MTRACE ("Instantiated.");
		}
		virtual ~calculate_timestep () {}
		
		virtual void execute () {
			plan::execute ();
			if (timestep != previous_timestep) {
				*flags_ptr &= ~factorized;
			}
			previous_timestep = timestep;
		};

	protected:
		double& timestep;
		double previous_timestep;
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
		bases::calculate_timestep::execute ();
	}

private:
	double initial_timestep;
};