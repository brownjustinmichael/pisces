/*!**********************************************************************
 * \file thermo_compositional.hpp
 * /Users/justinbrown/Dropbox/pisces/src/elements/data/thermo_compositional.hpp
 ************************************************************************/

#include "implemented_data.hpp"

namespace data
{
	/*!**********************************************************************
	 * \brief A data object designed to hold and output thermo-compositional data
	 ************************************************************************/
	class thermo_compositional_data : public implemented_data
	{
	protected:
		using implemented_data::n;
		using implemented_data::m;
		using implemented_data::grid_m;
		using implemented_data::iterator;
		using implemented_data::duration;
		
		std::vector <double> area; //!< A vector containing the area of each cell, for weighted averages
		
	public:
		using implemented_data::initialize;

		/*!**********************************************************************
		 * \param i_axis_n The horizontal axis object
		 * \param i_axis_m The vertical axis object
		 * \param i_name The integer name of the element
		 * \param n_elements The total number of elements
		 * \param i_params The parameters object associated with the run
		 ************************************************************************/
		thermo_compositional_data (grids::axis *i_axis_n, grids::axis *i_axis_m, int i_name, int n_elements, io::parameters& i_params);
		
		virtual ~thermo_compositional_data () {}

		static std::shared_ptr<data> instance(grids::axis* i_axis_n, grids::axis* i_axis_m, int i_name, int n_elements, io::parameters& i_params);
	};
} /* data */