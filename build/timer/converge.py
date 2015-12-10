from timing import Timer, Argument

mpiProcs = Argument ("mpirun -np %d", extent = [1], prepend = True, processes = True)
maxThreads = Argument ("-V parallel.maxthreads %d", extent = [8], threads = True)
points = Argument ("-V grid.z.points %d", extent=[512,256,128,64,32])

timer = Timer ("isces", 
               mpiProcs,
               maxThreads, 
               points,
               Argument ("-V grid.x.points %d", value=1) * points,
               setupCommandRoot = "isces_init", 
               directory = "../../run", 
               commandArgs = ["-D2", "-V", "output.stat.output", "true", "-V", "output.stat.every", "1"],
               uniques = [Argument ("-V input.file input_%03d_%%02i"), Argument ("-V output.stat.file stat_%03d_%%02i", name = "stat")])

results = timer.calculateTimes ()
