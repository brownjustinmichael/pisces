from timing import Timer, Argument

mpiProcs = Argument ("mpirun -np %d", extent = [1, 2, 4, 8], prepend = True, processes = True)
maxThreads = Argument ("-V parallel.maxthreads %d", extent = [1, 2, 4, 8], threads = True)

timer = Timer ("pisces", 
               mpiProcs,
               maxThreads, 
               Argument ("-V parallel.transform.threads %d", extent = [1, 2, 4, 8], upperBound = maxThreads), 
               setupCommandRoot = "pisces_init", 
               directory = "../../run", 
               commandArgs = ["-D2"], 
               uniques = [Argument ("-V input.file input_%03d_%%02i")])

results = timer.calculateTimes (torque = False, iterations = 2)

print (results)