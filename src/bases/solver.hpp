/*!***********************************************************************
 * \file bases/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_BDH126SH
#define SOLVER_HPP_BDH126SH
	
#include "plan.hpp"
#include "transform.hpp"

/*!*******************************************************************
 * \brief Execution flags used by the solver class
 *********************************************************************/
enum solver_flags {
	factorized = 0x08,
	never_factorized = 0x04
};

namespace bases
{
	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 *********************************************************************/
	class solver : public explicit_plan
	{
	public:
		/*!*******************************************************************
		 * \param i_flags_ptr A pointer to the integer execution flags
		 *********************************************************************/
		solver (int i_n, int i_name_in, int i_name_out) : explicit_plan (i_n, i_name_in, i_name_out) {}
		
		virtual ~solver () {}
			
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 *********************************************************************/
		virtual void execute () {
			explicit_plan::execute ();
			if ((!(*flags_ptr & factorized)) || (*flags_ptr & never_factorized)) {
				factorize ();
			}
			*flags_ptr |= transformed;
		}
		
	protected:		
		/*!*******************************************************************
		 * \brief Factorize the matrix equation
		 * 
		 * This method checks first whether the matrix has been factorized, 
		 * according to the execution flags. Upon success, it notes in the
		 * execution flags that the matrix has been factorized.
		 *********************************************************************/
		virtual void factorize () {
			*flags_ptr |= factorized;
		}
	};
} /* bases */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
