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

dims=tuple()

i=0
for dimension in data.dimensions:
	if i!=axis:
		new.createDimension(dimension, len(data.dimensions[dimension]))
		dims+=(dimension,)
	i+=1

for variable in data.variables:
	new.createVariable(variable, data[variable].dtype, dimensions=dims)
	new[variable][:]=np.mean(data[variable], axis=axis)

# data.close()
new.close()