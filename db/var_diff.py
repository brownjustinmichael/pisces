import pisces_db as db
import matplotlib.pyplot as plt
import numpy

session=db.Session()

base=db.SimulationEntry.from_params("config_var.yaml")
q=session.query(db.SimulationEntry.Table)

var="flux_temperature"

ain = 1.0
prcz = 0.1
racz = 1.e7

for prrz in (0.1, 0.01, 0.001):
	sim=base.clone()
	sim.entry.grid__x__points=2**np
	sim.entry.grid__z__points=2**np
	# sim.entry.time__max=10.**dt
	# sim.entry.time__init=10.**dt
	sim.entry.time__cfl=float(10.**dt)
	sim.entry.output__stat__file="stat_prrz=%f_%%02i_%%%%02i" % (prrz)
	sim.entry.output__cart__file="cart_prrz=%i_%i_%%02i_%%%%02i" % (prrz)

	sim.entry.equations__z_velocity__sources__temperature = racz / prcz
	sim.entry.equations__temperature__korre_Ts__Ain = ain * prcz / prrz
	sim.entry.equations__temperature__korre_diff__rz_diffusion = 1. / prrz

	query=sim.same_sub(q, "equations")
	query=sim.same_sub(query, "time")
	query=sim.same_sub(query, "grid")

	if query.first () is None:
		session.add(sim.entry)
		print ("Running with Pr_RZ = " % (prrz))
		new=sim.run(cwd="../sims/variable", execute="../../run/isces", reset=False, debug=2)
		session.commit()

	print (query.first ().id)

	step=db.SimulationEntry(query.first()).steps(session).order_by(db.StepEntry.Table.step.desc()).first()
	print (step)
	dts.append(dt)
	nps.append(np)
	vals.append(getattr(step, var))

vals=numpy.abs ((numpy.array(vals) - vals[0])/vals[0])
dx=base.entry.grid__z__width / 2 ** numpy.array(nps)
dx=dx[vals>0.]
dts=2.**numpy.array(dts)[vals>0.]
vals=vals[vals>0.]
s=plt.scatter(dts, vals, c=dx)

print(vals)
print(dts)

plt.xscale("log")
plt.yscale("log")
plt.ylim((10**numpy.floor(numpy.log10(min(vals))), 10**numpy.ceil(numpy.log10(max(vals)))))
plt.xlim((10**numpy.floor(numpy.log10(min(dts))), 10**numpy.ceil(numpy.log10(max(dts)))))

plt.colorbar(s)

plt.show()

