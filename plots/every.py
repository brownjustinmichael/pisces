import netCDF4 as nc
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--every", default=10)
parser.add_argument("--axis", default="time")

args=parser.parse_args()

data=nc.Dataset(args.input)
new=nc.Dataset(args.output, "w")

for dimension in data.dimensions:
	if dimension == args.axis:
		new.createDimension(dimension, len(data.dimensions[dimension]) // int(args.every) + (1 if len(data.dimensions[dimension]) % int(args.every) != 0 else 0))
	else:
		new.createDimension(dimension, len(data.dimensions[dimension]))

for variable in data.variables:
	print("Extracting %s" % variable)
	var_dims=list(data[variable].dimensions)

	new.createVariable(variable, data[variable].dtype, dimensions=var_dims)

	new[variable][:]=data[variable][::int(args.every)]

# assert(False)
# data.close()
new.close()