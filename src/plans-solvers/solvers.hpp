/*!**********************************************************************
 * \file plans-solvers/solvers.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVERS_H__
#define SOLVERS_H__

#include "boundaries/implemented_boundary.hpp"
#include "equation.hpp"
#include "solvers/collocation.hpp"
#include "solvers/fourier.hpp"
#include "solvers/incompressible.hpp"
#include "solvers/pseudo.hpp"

namespace plans
{
	/**
	 * @brief Shorthand to generate a split directional solve
	 * 
	 * @param equation_ptr A shared pointer to an equation
	 * @param timestep A reference to the timestep duration of the system
	 * @param boundary_0 A boundary factory to the top boundary
	 * @param boundary_n A boundary factory to the bottom boundary
	 * @return The equation pointer for use in equation definitions
	 */
	template <class datatype>
	std::shared_ptr <solvers::equation <datatype>> split_solver (std::shared_ptr <typename solvers::equation <datatype>> equation_ptr, datatype &timestep, std::shared_ptr <typename boundaries::boundary <datatype>::factory> boundary_0, std::shared_ptr <typename boundaries::boundary <datatype>::factory> boundary_n) {
		if (equation_ptr->messenger_ptr->get_id () > 0) {
			boundary_0 = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::communicating_boundary <datatype>::factory (equation_ptr->messenger_ptr));
		}
		if (equation_ptr->messenger_ptr->get_id () + 1 < equation_ptr->messenger_ptr->get_np ()) {
			boundary_n = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::communicating_boundary <datatype>::factory (equation_ptr->messenger_ptr));
		}
		equation_ptr->add_solver (typename solvers::collocation <datatype>::factory (equation_ptr->messenger_ptr, timestep, *boundary_0, *boundary_n), solvers::z_solver);
		equation_ptr->add_solver (typename solvers::fourier <datatype>::factory (timestep, *boundary_0, *boundary_n), solvers::x_solver);
		return equation_ptr;
	}

	/**
	 * @brief Shorthand to generate a divergence corrector
	 * 
	 * @param equation_ptr A shared pointer to an equation
	 * @param vel_n_ptr A pointer to the x directional velocity equation
	 * @param vel_m_ptr A pointer to the z directional velocity equation
	 * @return The equation pointer for use in equation definitions
	 */
	template <class datatype>
	std::shared_ptr <solvers::equation <datatype>> div (std::shared_ptr <solvers::equation <datatype>> equation_ptr, std::shared_ptr <solvers::equation <datatype>> vel_n_ptr, std::shared_ptr <solvers::equation <datatype>> vel_m_ptr) {
		std::shared_ptr <typename boundaries::boundary <datatype>::factory> boundary_0, boundary_n;
		if (equation_ptr->messenger_ptr->get_id () > 0) {
			boundary_0 = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::communicating_boundary <datatype>::factory (equation_ptr->messenger_ptr));
		}
		if (equation_ptr->messenger_ptr->get_id () + 1 < equation_ptr->messenger_ptr->get_np ()) {
			boundary_n = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::communicating_boundary <datatype>::factory (equation_ptr->messenger_ptr));
		}
		equation_ptr->add_solver (typename solvers::incompressible <datatype>::factory (equation_ptr->messenger_ptr,  *boundary_0, *boundary_n, *vel_n_ptr, *vel_m_ptr), solvers::x_solver);
			return equation_ptr;
	}

	/**
	 * @brief Shorthand to generate a divergence corrector
	 * 
	 * @param equation_ptr A shared pointer to an equation
	 * @param vel_n_ptr A pointer to the x directional velocity equation
	 * @param vel_m_ptr A pointer to the z directional velocity equation
	 * @param density A pointer to the density data
	 * @param pressure A pointer to the pressure data
	 * @param gamma The ratio of specific heats
	 * @return The equation pointer for use in equation definitions
	 */
	template <class datatype>
	std::shared_ptr <solvers::equation <datatype>> pdiv (std::shared_ptr <solvers::equation <datatype>> equation_ptr, std::shared_ptr <solvers::equation <datatype>> vel_n_ptr, std::shared_ptr <solvers::equation <datatype>> vel_m_ptr, datatype *density, datatype *pressure, datatype gamma = 5./3.) {
		std::shared_ptr <typename boundaries::boundary <datatype>::factory> boundary_0, boundary_n;
		if (equation_ptr->messenger_ptr->get_id () > 0) {
			boundary_0 = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::communicating_boundary <datatype>::factory (equation_ptr->messenger_ptr));
		}
		if (equation_ptr->messenger_ptr->get_id () + 1 < equation_ptr->messenger_ptr->get_np ()) {
			boundary_n = std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::communicating_boundary <datatype>::factory (equation_ptr->messenger_ptr));
		}
		equation_ptr->add_solver (typename solvers::pseudo_incompressible <datatype>::factory (equation_ptr->messenger_ptr,  *boundary_0, *boundary_n, *vel_n_ptr, *vel_m_ptr, density, pressure, gamma), solvers::x_solver);
			return equation_ptr;
	}

	/**
	 * @brief Shorthand to generate dirichlet boundary conditions
	 * 
	 * @param value The value to be set at the dirichlet boundary
	 * @return The boundary factory associated with the dirichlet boundary
	 */
	template <class datatype>
	typename std::shared_ptr <typename boundaries::boundary <datatype>::factory> dirichlet (datatype value) {
		return std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::fixed_boundary <datatype>::factory (value));
	}

	/**
	 * @brief Shorthand to generate neumann boundary conditions
	 * 
	 * @param value The value to be set at the neumann boundary
	 * @return The boundary factory associated with the neumann boundary
	 */
	template <class datatype>
	typename std::shared_ptr <typename boundaries::boundary <datatype>::factory> neumann (datatype value) {
		return std::shared_ptr <typename boundaries::boundary <datatype>::factory> (new typename boundaries::fixed_deriv_boundary <datatype>::factory (value));
	}
}

#endif /* end of include guard: SOLVERS_H__ */
