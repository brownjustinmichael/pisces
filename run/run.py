import os
import argparse

import pisces_utils.config as config
import pisces_utils.launch as launch
import pisces_utils.torque_launch as tlaunch

parser = argparse.ArgumentParser()

parser.add_argument("config_file", nargs="?", default="config.yaml")
parser.add_argument("--np", type=int, default=1)
parser.add_argument("--init", default=None)
parser.add_argument("--code", default="ISCES")
parser.add_argument("--launcher", default="Launcher")
parser.add_argument("--mpi", default="mpirun")

args = parser.parse_args()

configuration = config.Configuration(args.config_file)
configuration["np"] = args.np

code = launch.CodeRegistry.registry[args.code](configuration, init=args.init)

launcher = launch.LauncherRegistry.registry[args.launcher](code)
launcher.launch(init=(args.init not in ["False", "false", "f", "F"]), mpi=args.mpi)
launcher.wait()