import netCDF4 as nc
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--axis", default="z")
parser.add_argument("--deriv", default=None)

args=parser.parse_args()

data=nc.Dataset(args.input)
new=nc.Dataset(args.output, "w")

axis=list(data.dimensions.keys()).index(args.axis)
deriv=args.deriv

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
		axis=var_dims.index(args.axis)
		var_dims.pop(axis)
		to_mean=True
	except ValueError:
		pass

	new.createVariable(variable, data[variable].dtype, dimensions=var_dims)

	if to_mean:
		new[variable][:]=np.mean(data[variable], axis=axis)
	else:
		new[variable][:]=data[variable][:]

	if deriv is not None and deriv in list(data[variable].dimensions):
		new.createVariable(variable + "_deriv", data[variable].dtype, dimensions=var_dims)
		new[variable + "_deriv"][:]=0.
		new[variable + "_deriv"][...,:-1]=np.mean(np.diff(data[variable], axis=list(data[variable].dimensions).index(deriv)), axis=axis)

# assert(False)
# data.close()
new.close()