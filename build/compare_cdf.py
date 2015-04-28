#! /usr/bin/env python

import netCDF4
import argparse
import numpy as np

parser = argparse.ArgumentParser ()
parser.add_argument ('file1')
parser.add_argument ('file2')
parser.add_argument ('--diff', default = 1.0e-8)

namespace = parser.parse_args ()

file1 = netCDF4.Dataset (namespace.file1)
file2 = netCDF4.Dataset (namespace.file2)

for key in file1.variables:
    diff = np.array (file2.variables [key]) - np.array (file1.variables [key])
    if np.max (np.abs (diff)) > namespace.diff:
        raise ValueError ("Files differ")
        
print ("Success!")
