from timing import Timer, Argument

mpiProcs = Argument ("mpirun -np %d", extent = [1,2], prepend = True, threaded = True)
maxThreads = Argument ("-V parallel.maxthreads %d", extent = [1, 2], threaded = True)

timer = Timer ("pisces", 
               mpiProcs,
               maxThreads, 
               Argument ("-V parallel.transform.threads %d", extent = [1,2], upperBound = maxThreads), 
               setupCommandRoot = "pisces_init", 
               directory = "../../run", 
               commandArgs = ["-D2"], 
               uniques = [Argument ("-V input.file input_%03d_%%02i")])

results = timer.calculateTimes (torque = True, iterations = 2)

print (results)