import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata, LinearNDInterpolator

class Output(nc.Dataset):
	"""An output object for analysis"""
	def __init__(self, file, *args, **kwargs):
		super(Output, self).__init__(file, *args, **kwargs)

	@classmethod
	def from_simulation_entry(cls, simulation_entry, **kwargs):
		files = simulation_entry.get_files(**kwargs)
		outputs = []
		for file in files:
			outputs.append(cls(file))
		return outputs

	def norm(self, other, variable, self_index=-1, other_index=-1, n=2, positions=("x", "z")):
		"""
		Calculate the Ln norm of the file with another. If the two are not the same size, self will be interpolated via nearest neighbors to the points of other.

		:type other: :class:`Output`
		:param other: The other output object from which the difference should be calculated
		:type variable: :class:`str`
		:param variable: The variable key to use for the norm calculation
		:type index: :class:`int`
		:param index: The index to use for the norm calculation
		:type n: :class:`int`
		:param n: The order of the norm
		:type positions: :class:`tuple` of :class:`str`
		:param positions: A tuples of variable keys to use as positional arguments to interpolate if necessary

		:return: The relative Ln norm difference between the two outputs
		"""
		try:
			diff = (self.variables[variable][self_index] - other.variables[variable][other_index])
		except ValueError:
			flattened = [self.variables[var][self_index].flatten() for var in positions]
			points = [i for i in zip(*flattened)]
			values = self.variables[variable][self_index].flatten()

			inter_flattened = [other.variables[var][other_index].flatten() for var in positions]
			inter_points = [i for i in zip(*inter_flattened)]

			interpolator = LinearNDInterpolator(points, values)
			new_values = np.array([interpolator(point) for point in inter_points])
			diff = (new_values - other.variables[variable][other_index].flatten())

		scale = np.sum(other.variables[variable][other_index])
		return (np.sum(diff ** n)) ** (1.0 / n) / scale
