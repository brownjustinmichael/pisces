from os.path import join
from datetime import datetime
from subprocess import call
import copy

class Timer (object):
    """
    An object that times the execution of a command given a set of possible variances
    """
    def __init__(self, commandRoot, *args, **kwargs):
        super(Timer, self).__init__()
        self.commandRoot = commandRoot
        self.setupCommandRoot = kwargs.pop ("setupCommandRoot", None)
        self.commandArgs = kwargs.pop ("commandArgs", [])
        self.variances = []
        self.directory = kwargs.pop ("directory", None)
        for arg in args:
            if not isinstance (arg, Variance):
                arg = Variance (*arg)
            self.variances.append (arg)
        self.current = self.variances [:]
            
    def getSetupCommand (self, *args, **kwargs):
        command = [join (self.directory, self.setupCommandRoot)] + self.commandArgs
        for arg in args:
            arg (command)
        for var in self.variances:
            if var not in args:
                var (command)
        # print (command)
        return command
        
    def getCommand (self, *args, **kwargs):
        command = [join (self.directory, self.commandRoot)] + self.commandArgs
        for arg in args:
            arg (command)
        for var in self.variances:
            if var not in args:
                var (command)
        # print (command)
        return command
                
    def time (self, *args):
        if self.setupCommandRoot is not None:
            call (self.getSetupCommand (*args))
        startTime = datetime.now()
        call (self.getCommand (*args))
        return datetime.now() - startTime
        
    def calculateMinTime (self, *args):
        baseTime = self.time ()
        
        times = {tuple (self.variances) : baseTime}
        currentTime = baseTime
        if len (self.variances) == 0:
            return
        tests = [(var,) for var in self.variances [0].generate ()]
        for variance in self.variances [1:]:
            tests = [i + (var,) for i in tests for var in variance.generate ()]
        for variances in tests:
            if variances not in times:
                currentTime = self.time (*variances)
                times [variances] = currentTime
        best = min(times, key=times.get)
        print (best, baseTime / times [best])
        return times

class Variance (object):
    """
    An object that can return a part of a command associated with varying a particular parameter
    """
    def __init__(self, command, extent = None, value = 1, lowerBound = 1, **kwargs):
        super(Variance, self).__init__()
        self.command = command
        self.extent = extent
        if value < lowerBound:
            raise RuntimeError ("Can't have a value below %s" % str (lowerBound))
        self.value = value
        self.lowerBound = lowerBound
        self.kwargs = kwargs
        
    def copy (self, **kwargs):
        newkwargs = copy.copy (self.kwargs)
        newkwargs ["value"] = self.value
        newkwargs ["lowerBound"] = self.lowerBound
        newkwargs ["extent"] = self.extent
        for arg in kwargs:
            newkwargs [arg] = kwargs [arg]
        newObject = Variance (self.command, **newkwargs)
        return newObject
        
    def __mul__ (self, scalar):
        return Variance (self.command, self.value * scalar, self.lowerBound, **self.kwargs)
        
    def __add__ (self, scalar):
        return Variance (self.command, self.value + scalar, self.lowerBound, **self.kwargs)
        
    def __div__ (self, scalar):
        return Variance (self.command, self.value / scalar, self.lowerBound, **self.kwargs)
    
    def __sub__ (self, scalar):
        return Variance (self.command, self.value - scalar, self.lowerBound, **self.kwargs)
        
    def __call__ (self, fullCommand = None):
        if fullCommand is None:
            fullCommand = []
        if self.kwargs.get ("prepend", False):
            for subcommand in (self.command + str (self.value)).split () [::-1]:
                fullCommand.insert (0, subcommand)
        else:
            for subcommand in (self.command + str (self.value)).split ():
                fullCommand.append (subcommand)
        return fullCommand
        
    def __eq__ (self, other):
        return self.command == other.command
        
    def __hash__ (self):
        return hash (self.command) + hash (self.value)
        
    def __str__ (self):
        return " ".join (self ())
    
    def __repr__ (self):
        return "<%s>" % str (self)
        
    def generate (self):
        if self.extent is None:
            yield self.copy ()
        else:
            for value in self.extent:
                yield self.copy (value = value)
        
timer = Timer ("pisces", Variance ("mpiexec -np ", extent = [1,2], prepend = True), Variance ("-V parallel.maxthreads ", extent = [1]), Variance ("-V parallel.transform.threads ", extent = [1,2]), setupCommandRoot = "pisces_init", directory = "../../run", commandArgs = ["-D7"])

results = timer.calculateMinTime ()