import netCDF4 as nc
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")

args=parser.parse_args()

data=nc.Dataset(args.input)
new=nc.Dataset(args.output, "w")

axis=1
deriv="z"

dims=[]

i=0
for dimension in data.dimensions:
	new.createDimension(dimension, len(data.dimensions[dimension]))
	if i!=axis:
		dims+=[dimension,]
	else:
		removed=dimension
	i+=1

for variable in data.variables:
	print("Extracting %s" % variable)
	var_dims=list(data[variable].dimensions)

	to_mean=False
	try:
		var_dims.pop(var_dims.index(removed))
		to_mean=True
	except ValueError:
		pass
	new.createVariable(variable, data[variable].dtype, dimensions=var_dims)
	if to_mean:
		new[variable][:]=np.mean(data[variable], axis=axis)
	else:
		new[variable][:]=data[variable][:]

	if deriv in var_dims:
		new.createVariable(variable + "_deriv", data[variable].dtype, dimensions=var_dims)
		new[variable + "_deriv"][:]=0.
		new[variable + "_deriv"][...,:-1]=np.diff(new[variable], axis=var_dims.index(deriv))

# assert(False)
# data.close()
new.close()