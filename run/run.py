#! /usr/bin/env python

import os
import argparse
import yaml

import pisces_utils.db as db
import pisces_utils.config as config
import pisces_utils.launch as launch

exports = {}
custom_run_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "run_custom.yaml")
if os.path.isfile(custom_run_file):
	new_exports = yaml.load(open(custom_run_file))
	if new_exports is not None:
		for key in new_exports:
		    exports [key] = new_exports [key]

parser = argparse.ArgumentParser()

parser.add_argument("config_file", nargs="?", default="config.yaml")
parser.add_argument("--np", type=int, default=1)
parser.add_argument("--init", default=None)
parser.add_argument("--from_dump", default=False, type=bool)
parser.add_argument("--code", default="PISCES")
parser.add_argument("--launcher", default="Launcher")
parser.add_argument("--debug", default=2, type=int)

args = parser.parse_args()

configuration = config.Configuration(args.config_file)
if args.np != 1 or "np" not in configuration:
	configuration["np"] = args.np

Code = launch.CodeRegistry.registry[args.code]
Launcher = launch.LauncherRegistry.registry[args.launcher]

code = Code(configuration, init=args.init)

session = db.Session()
entry = db.SimulationEntry.query(session, **configuration)

if entry is not None and entry.date >= code.date:
	print("Up-to-date db entry already exists for simulation.")
	code = Code.from_simulation_entry(entry)
	args.from_dump = True
session.close()

launcher = Launcher(code)
try:
	print("Launching simulation")
	launcher.launch(from_dump=args.from_dump, **exports)
	launcher.wait()
except KeyboardInterrupt:
	print("Received keyboard interrupt. Canceling task.")
	launcher.cancel()

print("Recording in database.")
session = db.Session()
code.record(session)
