#! /usr/bin/env python

import os
import argparse
import yaml

import pisces_utils.config as config
import pisces_utils.launch as launch

exports = {}
try:
	with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "run_custom.yaml"), 'r') as stream:
		new_exports = yaml.load(stream)
		for key in new_exports:
		    exports [key] = new_exports [key]
except IOError:
	pass
except TypeError:
	pass

parser = argparse.ArgumentParser()

parser.add_argument("config_file", nargs="?", default="config.yaml")
parser.add_argument("--code", default="PISCES")

args = parser.parse_args()

configuration = config.Configuration(args.config_file)

code = launch.CodeRegistry.registry[args.code](configuration)

code.record()