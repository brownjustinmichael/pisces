import pisces_db as db
import matplotlib.pyplot as plt
import numpy

session=db.Session()

base=db.SimulationEntry.from_params("config_orig.yaml")
q=session.query(db.SimulationEntry.Table)
q=base.same_sub(q, "equations")

var="flux_temperature"
dts=[]
nps=[]
vals=[]

for dt in range(2, -5, -1):
	for np in range(6, 11):
		sim=base.clone()
		sim.entry.grid__x__points=2**np
		sim.entry.grid__z__points=2**np
		sim.entry.time__max=2.**dt
		sim.entry.time__init=2.**dt
		sim.entry.output__stat__file="stat_%i_%i_%%02i_%%%%02i" % (dt, np)
		query=sim.same_sub(q, "time")
		query=sim.same_sub(query, "grid")

		if query.first() is None or db.SimulationEntry(query.first()).steps(session).count() == 0:
			session.add(sim.entry)
			session.commit()
			new=sim.run(cwd="../sims/pisces_test/", session=session)
			session.commit()

		step=db.SimulationEntry(query.first()).steps(session).order_by(db.StepEntry.Table.step.desc()).first()
		dts.append(dt)
		nps.append(np)
		vals.append(getattr(step, var))

vals=numpy.abs ((numpy.array(vals) - vals[-1])/vals[-1])
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

