/*!***********************************************************************
 * \file bases/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_BDH126SH
#define SOLVER_HPP_BDH126SH

#include "transform.hpp"
	
/*!*******************************************************************
 * \brief Execution flags used by the solver class
 *********************************************************************/
enum solver_flags {
	factorized = 0x02,
	never_factorized = 0x04
};

namespace bases
{
	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 *********************************************************************/
	class solver
	{
	public:
		/*!*******************************************************************
		 * \param i_flags_ptr A pointer to the integer execution flags
		 *********************************************************************/
		solver (int *i_flags_ptr = NULL, int i_logger = -1) {
			logger = i_logger;
			if (i_flags_ptr) {
				flags_ptr = i_flags_ptr;
			} else {
				default_flags = never_factorized;
				flags_ptr = &default_flags;
			}
		}
		
		virtual ~solver () {}
			
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 *********************************************************************/
		virtual void solve () {
			factorize ();
			i_solve ();
			*flags_ptr |= transformed;
		}
		
	protected:
		int logger;
		int *flags_ptr; //!< A pointer to the integer execution flags
		int default_flags; //!< Default flags to use if i_flags_ptr is NULL
		
		/*!*******************************************************************
		 * \brief The implementation for factorizing the matrix equation
		 *********************************************************************/
		virtual void i_factorize () = 0;
	
		/*!*******************************************************************
		 * \brief The implementation for solving a factorized matrix equation
		 *********************************************************************/
		virtual void i_solve () = 0;
		
	private:
		/*!*******************************************************************
		 * \brief Factorize the matrix equation
		 * 
		 * This method checks first whether the matrix has been factorized, 
		 * according to the execution flags. Upon success, it notes in the
		 * execution flags that the matrix has been factorized.
		 *********************************************************************/
		virtual void factorize () {
			if ((! (*flags_ptr & factorized)) || (*flags_ptr & never_factorized)) {
				i_factorize ();
				*flags_ptr |= factorized;
			}
		}
	};
} /* bases */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
