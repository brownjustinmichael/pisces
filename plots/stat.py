#! /usr/bin/env python

import netCDF4 as nc
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--vel", default="z_velocity")

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

x_var = data.variables["x"]
z_var = data.variables["z"]
dx = np.diff(x_var, axis=list(data.dimensions.keys()).index("x")).resize(x_var.shape)
dz = np.diff(z_var, axis=list(data.dimensions.keys()).index("z")).resize(z_var.shape)
weight = dx * dz / np.sum(dx * dz)

for variable in data.variables:
	var_dims = list(data[variable].dimensions)

	mean_tuple = [i for i in range(1, len(var_dims))]

	print("Extracting %s" % variable)
	new.createVariable("max_" + variable, data[variable].dtype, dimensions=data.dimensions[0:])
	new.createVariable("avg_" + variable, data[variable].dtype, dimensions=data.dimensions[0:])
	new.createVariable("flux_" + variable, data[variable].dtype, dimensions=data.dimensions[0:])

	new["max_" + variable][:] = np.max(data[variable], axis=mean_tuple)
	new["avg_" + variable][:] = np.mean(data[variable] * weight, axis=mean_tuple)
	new["flux_" + variable][:] = np.mean(data[variable] * weight * data[args.vel], axis=mean_tuple)

new.close()