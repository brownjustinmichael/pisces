import pisces_db as db
import matplotlib.pyplot as plt
import numpy

session=db.Session()

base=db.SimulationEntry.from_params("config_convect.yaml")
q=session.query(db.SimulationEntry.Table)
q=base.same_sub(q, "equations")

var="flux_temperature"
dts=[]
nps=[]
vals=[]

dt=-4
np=8
for r in numpy.arange(0, 1, 0.1):
	sim=base.clone()
	sim.entry.grid__x__points=2**np
	sim.entry.grid__z__points=2**np
	sim.entry.time__max=2**dt
	sim.entry.time__init=2**dt
	print(sim.entry.equations__velocity__diffusion*657.5*10.**r)
	sim.entry.equations__z_velocity__sources__temperature=float(sim.entry.equations__velocity__diffusion*657.5*10.**r)
	sim.entry.output__stat__file="stat_%f_%%02i_%%%%02i" % (r)
	sim.entry.output__cart__file="cart_%f_%%02i_%%%%02i" % (r)
	query=sim.same_sub(q, "time")
	query=sim.same_sub(q, "grid")

	if query.first() is None:
		session.add(sim.entry)
		new=sim.run(cwd="../sims/convect/", execute="../../run/isces")
		session.commit()

	step=db.SimulationEntry(query.first()).steps.order_by(db.StepEntry.Table.step.desc()).first()
	dts.append(dt)
	nps.append(np)
	vals.append(getattr(step, var))

plt.scatter(2 ** numpy.array(nps), numpy.abs (numpy.array(vals) - vals[0]))

plt.xscale("log")
plt.yscale("log")

plt.show()

