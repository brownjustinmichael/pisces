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

/*!*******************************************************************
 * \brief Execution flags used by the solver class
 *********************************************************************/
enum solver_flags {
	factorized = 0x08,
};

namespace bases
{
	class element;

	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 *********************************************************************/
	class solver : public explicit_plan
	{
	public:
		/*!*******************************************************************
		 * \copydoc explicit_plan::explicit_plan ()
		 *********************************************************************/
		solver (element* i_element_ptr, int i_n, int i_name_in, int i_name_out, int i_flags = 0x00) : 
		explicit_plan (i_element_ptr, i_n, i_name_in, i_name_out, i_flags) {}
		
		virtual ~solver () {}
			
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 *********************************************************************/
		virtual void execute () {
			explicit_plan::execute ();
			if ((!(flags & factorized)) || !(*flags_ptr & factorized)) {
				factorize ();
			}
			
			/*
				TODO The transformed flag requires better handling. It could be handled entirely by element?
			*/
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
			flags |= factorized;
		}
	};
} /* bases */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
