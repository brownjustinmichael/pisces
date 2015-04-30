from timing import Timer, Argument

mpiProcs = Argument ("mpirun -np %d", extent = [1, 2, 4, 8], prepend = True, processes = True, runOnly = "-machinefile $HOSTFILE")
maxThreads = Argument ("-V parallel.maxthreads %d", extent = [8], threads = True)

timer = Timer ("pisces", 
               mpiProcs,
               maxThreads, 
               Argument ("-V parallel.transform.threads %d", extent = [1, 2], upperBound = maxThreads),
               Argument ("-V grid.z.points %d", value = 256) * mpiProcs, 
               setupCommandRoot = "pisces_init", 
               directory = "../../run", 
               commandArgs = ["-D2"],
               uniques = [Argument ("-V input.file input_%03d_%%02i")])

# results = timer.calculateTimes (torque = True, iterations = 4, hours = 4)

# print (results)