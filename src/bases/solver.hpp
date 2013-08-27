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
	first_run = 0x10
};

namespace bases
{
	/*!*******************************************************************
	 * \brief A class designed to solve a matrix equation
	 *********************************************************************/
	template <class datatype>
	class solver : public explicit_plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \copydoc explicit_plan::explicit_plan ()
		 *********************************************************************/
		solver (element <datatype>* i_element_ptr, int i_n, datatype* i_data_in, datatype* i_data_out, int i_flags = 0x00) : 
		explicit_plan <datatype> (i_element_ptr, i_n, i_data_in, i_data_out, i_flags) {}
		
		virtual ~solver () {}
			
		/*!*******************************************************************
		 * \brief Solve the matrix equation
		 *********************************************************************/
		virtual void execute () {
			if ((!(flags & factorized)) || !(*flags_ptr & factorized)) {
				_factorize ();
			}
			explicit_plan <datatype>::execute ();
			
			/*
				TODO The transformed flag requires better handling. It could be handled entirely by element?
			*/
		}
		
		virtual void factorize () {
			if ((!(flags & factorized)) || !(*flags_ptr & factorized)) {
				_factorize ();
			}
			flags |= factorized;
		}
		
	protected:
		/*!*******************************************************************
		 * \brief Factorize the matrix equation
		 * 
		 * This method checks first whether the matrix has been factorized, 
		 * according to the execution flags. Upon success, it notes in the
		 * execution flags that the matrix has been factorized.
		 *********************************************************************/
		virtual void _factorize () = 0;
		
		using explicit_plan <datatype>::flags;
		using explicit_plan <datatype>::flags_ptr;
	};
} /* bases */

#endif /* end of include guard: SOLVER_HPP_BDH126SH */
