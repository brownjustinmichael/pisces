import pisces_db as db
import matplotlib.pyplot as plt
import numpy

session=db.Session()

base=db.SimulationEntry.from_params("config_convect.yaml")
q=session.query(db.SimulationEntry.Table)

var="flux_temperature"
dts=[]
nps=[]
vals=[]

for dt in numpy.arange(-0.5, -2.0, -0.25):
	for np in range(4, 8):
		print ("Checking with dt = %f, np = %i" % (10.**dt, 2**np))
		sim=base.clone()
		sim.entry.grid__x__points=2**np
		sim.entry.grid__z__points=2**np
		# sim.entry.time__max=10.**dt
		# sim.entry.time__init=10.**dt
		sim.entry.time__cfl=float(10.**dt)
		sim.entry.output__stat__file="stat_low_%f_%i_%%02i_%%%%02i" % (dt, np)
		sim.entry.output__cart__file="cart_low_%f_%i_%%02i_%%%%02i" % (dt, np)
		query=sim.same_sub(q, "equations")
		query=sim.same_sub(query, "time")
		query=sim.same_sub(query, "grid")

		if query.first () is None:
			session.add(sim.entry)
			print ("Running with dt = %f, np = %i" % (10.**dt, 2**np))
			new=sim.run(cwd="../sims/convect/", execute="../../run/isces", reset=False, debug=3)
			session.commit()

		print (query.first ().id)

		step=db.SimulationEntry(query.first()).steps(session).order_by(db.StepEntry.Table.step.desc()).first()
		print (step)
		dts.append(10.**dt)
		nps.append(np)
		vals.append(getattr(step, var))

		scale_vals=numpy.abs ((numpy.array(vals) - vals[-1])/vals[-1])
		# vals=numpy.abs ((numpy.array(vals)))
		dx=base.entry.grid__z__width / 2 ** numpy.array(nps)
		scale_dx=dx[scale_vals>0.]
		scale_dts=numpy.array(dts)[scale_vals>0.]
		scale_vals=scale_vals[scale_vals>0.]
		
		if len (scale_vals) > 0:
			plt.clf ()

			plt.xscale("log")
			plt.yscale("log")

			# plt.ylim((10**numpy.floor(numpy.log10(min(scale_vals))), 10**numpy.ceil(numpy.log10(max(scale_vals)))))
			# plt.xlim((10**numpy.floor(numpy.log10(min(scale_dts))), 10**numpy.ceil(numpy.log10(max(scale_dts)))))

			s=plt.scatter(scale_dts, scale_vals, c=scale_dx)
			plt.colorbar(s)

			plt.draw ()
			plt.pause(0.0001) 


plt.show ()


# [  6.21590946e-01   5.72423735e-02   1.42262605e-05   2.14029603e-01   5.71841510e-02   3.10420927e-06   2.13842404e-01   5.71581537e-02]
# [ 0.70710678  0.70710678  0.70710678  0.5         0.5         0.5  0.35355339  0.35355339]



