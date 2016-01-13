import os
import argparse

import pisces_utils.config as config
import pisces_utils.launch as launch
import pisces_utils.torque_launch as tlaunch

parser = argparse.ArgumentParser()

parser.add_argument("config_file", nargs="?", default="config.yaml")
parser.add_argument("--np", dtype=int, default=1)
parser.add_argument("--init", default="init")
parser.add_argument("--lancher", default="Launcher")

args = parser.parse_args()

configuration = config.Configuration(args.config_file)
configuration.np = args.np

code = launch.ISCES(configuration, init=args.init)

launcher = LauncherRegistry.registry[args.launcher](code)
launcher.launch()
launcher.wait()