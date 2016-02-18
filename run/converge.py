import os

import pisces_utils.db as db
import pisces_utils.launch as launch
from pisces_utils.analysis.output import Output
from pisces_utils.config import Configuration
import matplotlib.pyplot as plt
import numpy
import yaml
import copy

session=db.Session()

dts=[]
nps=[]
vals=[]

base_yaml = Configuration("converge.yaml")
code = launch.PISCES(base_yaml)
launcher = launch.Launcher(code)
cd = os.getcwd()

def run(np, dt):
	os.chdir(cd)
	# print("Checking with dt = %f, np = %i" % (dt, np))
	code.config["grid"]["x"]["points"] = int(np)
	code.config["grid"]["z"]["points"] = int(np)
	code.config["time"]["max"] = float(dt)
	code.config["time"]["init"] = float(dt)

	code.config["wd"]=os.path.abspath(os.path.join(cd, "../sims/converge/diffusion/diff_%04d_%f" % (np, dt)))

	entry = db.SimulationEntry.query(session, **code.config)
	print(entry)
	if entry is None or entry.date is None or entry.date < code.date:
		print("Running with dt = %f, np = %i" % (dt, np))
		launcher.launch(check_dir=False)
		launcher.wait()
	entry = code.record(session)
	if entry is None:
		raise RuntimeError("Something weird happened")
	return entry

base = run(1024, 0.001)
solution = Output.from_simulation_entry(base)[-1]
runs = []
nps = []
dts = []
results = []
for dt in numpy.arange(2, 10, 2):
	for np in range(4, 10):
		entry = run(2**np, 1. / dt)
		runs.append(entry)
		output = Output.from_simulation_entry(entry)[-1]
		nps.append(2**np)
		dts.append(1. / dt)
		results.append(solution.norm(output, variable="scalar"))

fig, ax = plt.subplots(1, 1)
s = ax.scatter(nps, results, c=numpy.log10(dts))
ax.set_xscale("log")
ax.set_xlim((10,1.e3))
ax.set_yscale("log")
ax.set_ylim((1.e-7,1.e0))
cb = fig.colorbar(s)

plt.show()