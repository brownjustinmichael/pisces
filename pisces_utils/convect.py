import pisces_db as db
import matplotlib.pyplot as plt
import numpy

session=db.Session()

base=db.SimulationEntry.from_params("config_convect.yaml")
q=session.query(db.SimulationEntry.Table)

var="flux_temperature"
dts=[]
rs=[]
vals=[]
prs=[]

dt=-2.5
np=7
for r in numpy.arange(-1., 1.1, 0.1):
	sim=base.clone()
	sim.entry.grid__x__points=2**np
	sim.entry.grid__z__points=2**np
	sim.entry.time__max=10.**dt
	sim.entry.time__init=10.**-10
	sim.entry.equations__z_velocity__sources__temperature=float(sim.entry.equations__velocity__diffusion*657.5*10.**r)
	sim.entry.output__stat__file="stat_%f_%%02i_%%%%02i" % (r)
	sim.entry.output__cart__file="cart_%f_%%02i_%%%%02i" % (r)
	query=sim.same_sub(q, "equations")
	query=sim.same_sub(query, "time")
	query=sim.same_sub(query, "grid")

	if query.first() is None:
		session.add(sim.entry)
		print ("Running with RaPr = %f" % 10.**r)
		new=sim.run(cwd="../sims/convect/", execute="../../run/isces", reset=True, debug=3)
		session.commit()

	step=db.SimulationEntry(query.first()).steps(session).order_by(db.StepEntry.Table.step.desc()).first()
	dts.append(dt)
	rs.append(10.**r)
	prs.append(sim.entry.equations__velocity__diffusion)
	vals.append(getattr(step, var))

prs=numpy.array(prs)
vals=numpy.array(vals)

print(vals)
plt.scatter(rs, (vals + 1), label="ISCES")

rbdata = numpy.genfromtxt ("../documentation/data/moore_weiss_rb.dat")
plt.plot(rbdata[:,0], rbdata [:,1], label="Moore & Weiss (1973)")

plt.xscale("log")
plt.yscale("log")

plt.show()

