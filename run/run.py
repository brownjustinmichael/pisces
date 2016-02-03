#! /usr/bin/env python

import os
import argparse
import yaml

import pisces_utils.config as config
import pisces_utils.launch as launch
import pisces_utils.torque_launch as tlaunch

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
parser.add_argument("--np", type=int, default=1)
parser.add_argument("--init", default=None)
parser.add_argument("--code", default="ISCES")
parser.add_argument("--launcher", default="Launcher")
parser.add_argument("--debug", default=2, type=int)

args = parser.parse_args()

configuration = config.Configuration(args.config_file)
if args.np != 1 or "np" not in configuration:
	configuration["np"] = args.np

code = launch.CodeRegistry.registry[args.code](configuration, init=args.init)

exports ["init"] = (args.init not in ["False", "false", "f", "F"])

launcher = launch.LauncherRegistry.registry[args.launcher](code)
launcher.launch(**exports)
launcher.wait()