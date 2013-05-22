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

/*!*******************************************************************
 * \namespace bases
 *********************************************************************/
namespace bases
{
	/*!*******************************************************************
	 * \brief A plan that calculates the next timestep
	 *********************************************************************/
	class calculate_timestep : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_timestep The double reference to the timestep
		 *********************************************************************/
		calculate_timestep (double& i_timestep) : timestep (i_timestep) {
			previous_timestep = 0.0;
			MTRACE ("Instantiated.");
		}
		virtual ~calculate_timestep () {}
		
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () {
			plan::execute ();
			if (timestep != previous_timestep) {
				*flags_ptr &= ~factorized;
			}
			previous_timestep = timestep;
		};

	protected:
		double& timestep; //!< A double reference to the timestep
		double previous_timestep; //!< The double previous timestep
	};
} /* bases */

/*!*******************************************************************
 * \brief A plan that sets the timestep to a constant
 *********************************************************************/
class constant_timestep : public bases::calculate_timestep
{
public:
	/*!*******************************************************************
	 * \param i_initial_timestep The double timestep constant value
	 * \copydoc bases::calculate_timestep::calculate_timestep ()
	 *********************************************************************/
	constant_timestep (double i_initial_timestep, double& i_timestep) : calculate_timestep (i_timestep) {
		initial_timestep = i_initial_timestep;
	}
	
	virtual ~constant_timestep () {}
	
	/*!*******************************************************************
	 * \copydoc bases::calculate_timestep::execute ()
	 *********************************************************************/
	virtual void execute () {
		timestep = initial_timestep;
		bases::calculate_timestep::execute ();
	}

private:
	double initial_timestep; //!< The double timestep constant value
};