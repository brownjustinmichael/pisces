#!/usr/bin/env python
# the Scientific Python netCDF 3 interface
# http://dirac.cnrs-orleans.fr/ScientificPython/
from netCDF4 import Dataset
# the 'classic' version of the netCDF4 python interface
# http://code.google.com/p/netcdf4-python/
#from netCDF4_classic import Dataset
from numpy import arange, dtype, concatenate, delete # array module from http://numpy.scipy.org

from pyshell.config import DottedConfiguration

import sys
import os

def merge (input_file_rule, range1, range2, output_file_rule):
    # Get all the files to merge
    for i in range2:
        output = Dataset (output_file_rule % i, 'w')
        
        files = []
        for j in range1:
            string = input_file_rule % j
            files.append (Dataset (string % i, 'r'))
        
        # Get all the variables from the first file
        variables = []
        for var in files [0].variables:
            if (files [0].variables [var].shape != ()):
                variables.append (var);
        
        # Calculate the shape of the final array
        extent = 0
        for file in files:
            if file == files [0]:
                extent += file.variables [variables [0]].shape [-1]
            else:
                extent -= 2
                extent += file.variables [variables [0]].shape [-1] - 1
        shape = list (file.variables [variables [0]].shape)
        shape [-1] = extent
        shape = tuple (shape)
        
        # Set up the output
        data = {}
        dim_names = ["z", "x", "y"] [:len (shape)]
        dim_names.append (dim_names.pop (0))
        dim_names.reverse ()
        dim_names = tuple (dim_names)
        for i in range (len (shape)):
            output.createDimension(dim_names [i], shape [len (shape) - 1 - i])
        for var in variables:
            data [var] = output.createVariable (var, files [0].variables [var].dtype.char, dim_names)
        
        # Write to the files
        index = 0
        for file in files:
            if file != files [0]:
                index -= 2
            for var in variables:
                # Check that all these exist in all files
                if var not in file.variables:
                    print "%s not found in %s" % (var, str (file))
                    raise IndexError
                # Concatenate the data along the last axis
                else:
                    # Write the entire first file
                    if file == files [0]:
                        data [var] [:file.variables [var].shape [-1]] = file.variables [var] [...,:].transpose ()
                    else:
                        # Ignore the last two columns of the previous file and the first of the next one
                        data [var] [index:index + file.variables [var].shape [-1] - 1] = file.variables [var] [..., 1:].transpose ()
            if file == files [0]:
                index += file.variables [var].shape [-1]
            else:
                index += file.variables [var].shape [-1] - 1
                
            file.close ()
        output.close ()

def main ():
    config = DottedConfiguration.fromfile ("../input/config.yaml")
    
    if len (sys.argv) < 2:
        print "Proper use is './merge_cdf.py number_of_cdfs_to_merge'"
        exit ()
    
    merge ("../output/" + config ["output.file"] + ".cdf", range (int (sys.argv [1])), range (int (config ["time.steps"] / config ["output.every"]) + 1), "../output/" + config ["output.merge_file"] + ".cdf")
    
main ()
