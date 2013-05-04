/*!***********************************************************************
 * \file one_d/element.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_3SURDTOH
#define ELEMENT_HPP_3SURDTOH

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "../config.hpp"
#include "../bases/boundary.hpp"
#include "../bases/element.hpp"
#include "../bases/plan.hpp"
#include "../utils/utils.hpp"
	
namespace one_d
{
	/*!*******************************************************************
	 * A 1D implementation of the element base class. This provides the
	 * storage, indexing facilities, and failsafe_dump output. The plans should be added in a 
	 * further subclass.
	 * 
	 * \brief \copybrief bases::element
	 *********************************************************************/
	class element : public bases::element
	{
	public:	
		/*!*******************************************************************
		 * \param i_n The number of data elements in each scalar
		 * \copydoc bases::element::element ()
		 *********************************************************************/
		element (std::string i_name, int i_n, int i_flags) : bases::element (i_name, i_flags) {
			n = i_n;
			
			cell.resize (i_n);
			for (int i = 0; i < i_n; ++i) {
				cell [i] = i;
			}
			
			failsafe_dump = std::make_shared <io::simple_output> (io::simple_output ("dump_" + name + ".dat", n, logger));
			failsafe_dump->append (&cell [0]);
		}
		
		virtual ~element () {}
		
		/*!*******************************************************************
		 * \copydoc bases::element::operator[] ()
		 *********************************************************************/
		inline double& operator[] (int name) {
			if (scalars [name].size () == (unsigned int) 0) {
				scalars [name].resize (n, 0.0);
				failsafe_dump->append ((*this) [name]);
			}
			return scalars [name] [0];
		}
		
		virtual double& operator () (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
		
		inline void execute_boundaries () {
			bases::element::execute_boundaries ();
		}
		
		inline void update () {			
			bases::element::update ();
		}
		
		inline void explicit_reset () {
			bases::element::explicit_reset ();
			std::map <int, std::vector <double>>::iterator iter;
			for (iter = scalars.begin (); iter != scalars.end (); ++iter) {
				if (iter->first < 0) {
					utils::scale (n, 0.0, &(iter->second) [0]);					
				}
			}
		}
		
		inline int get_boundary_index (int edge) {
			MTRACE ("Getting boundary index...");
			if (edge == fixed_0) {
				return 0;
			} else if (edge == fixed_n) {
				return n - 1;
			} else {
				FATAL (logger, "Edge is not a one_d edge index.");
				throw 0;
			}
		}
		
		friend class boundary;
		
		bases::element::iterator begin () {return scalars.begin ();}
		
		bases::element::iterator end () {return scalars.end ();}
		
	protected:
		int n; //!< The number of elements in each 1D array
		std::vector<int> cell; //!< An integer array for tracking each cell number for output

		std::map <int, std::vector <double>> scalars; //!< A vector of scalar vectors
	};

	namespace chebyshev
	{
	/*!*******************************************************************
	 * \brief A simple implementation of the element class with diffusion
	 * 
	 * This class contains a full element's capacity to run a single 
	 * element diffusion in 1D. It cannot yet tie to other elements or 
	 * calculate its own timestep. 
	 *********************************************************************/
	class advection_diffusion_element : public element
	{
	public:
		/*!*******************************************************************
		 * \param i_n The number of elements in each 1D data array
		 * \param i_flags Flags for the boundary conditions and evaluation
		 *********************************************************************/
		advection_diffusion_element (std::string i_name, int i_n, double initial_position, double *intial_velocity, int i_flags,std::map<std::string,io::types>& inputParams);
		virtual ~advection_diffusion_element () {}
		
		inline void implicit_reset () {
			element::implicit_reset ();
			
			if (!(flags & factorized)) {
				utils::scale (n * n, 0.0, &matrix [0]);
			}
		}
		
		inline void update () {			
			element::update ();
		}
		
	private:
		std::vector<double> matrix; //!< A vector containing the double matrix used in the implicit solver
	};
	} /* chebyshev */
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
