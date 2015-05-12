from os.path import join
from datetime import datetime
from subprocess import Popen, call
import copy
from celery import Celery
import celeryconfig
import json
import operator

import socket

app = Celery('pisces_timer', backend = 'amqp')
app.config_from_object (celeryconfig)

@app.task
def timeCommand (command, setupCommand = None, iterations = 1, wrapperFile = "wrapper.py", processes = 1, threads = 1, torque = False, commandRoot = "job", hours = 1):
    if isinstance (command, str):
        command = [command]
    
    if setupCommand is not None:
        call (setupCommand)
        
    server_socket = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
    guess = 5000 + abs (hash (timeCommand.request.id)) % 1000
    
    while True:
        try:
            server_socket.bind (("", guess))
            break
        except socket.error:
            guess += 1
            
    server_socket.listen (5)
    
    print ("Sending to wrapper file...", wrapperFile)
    
    if torque:
        print ("Writing to batch file")

        batch_file = open ("batch_%04d.pbs" % guess, "w")

        batch_file.write ("#PBS -S /bin/bash\n")
        batch_file.write ("#PBS -q normal\n")
        batch_file.write ("#PBS -N %s\n" % commandRoot)

        batch_file.write ("#PBS -l nodes=%d:ppn=%d\n" % (processes, threads))
        batch_file.write ("#PBS -l walltime=%02i:00:00\n" % hours)
        batch_file.write ("cd $PBS_O_WORKDIR\n")
        batch_file.write ("cp $PBS_NODEFILE .\n")
        batch_file.write ("export HOSTFILE=hostfile_%04d\n" % guess)
        batch_file.write ("export I_MPI_PIN_DOMAIN=omp")
        batch_file.write ("export KMP_AFFINITY=compact")
        
        batch_file.write ("module load python\n")
        batch_file.write ("module switch python python/3.4.1\n")
        batch_file.write ("python3 host_rewrite.py $PBS_NODEFILE --ppn %d > $HOSTFILE\n" % (threads))
        # batch_file.write ("OMP_NUM_THREADS=%d\n" % threads)
        batch_file.write ("export OMP_NUM_THREADS=%d\n" % threads)

        batch_file.write (" ".join (["python", wrapperFile, "-a", str (socket.gethostbyname(socket.gethostname())), "-p", str (guess)]))
        batch_file.write ("\n")
        batch_file.close ()
        
        Popen (["qsub", "batch_%04d.pbs" % guess])
    else:
        Popen (["python3", wrapperFile, "-a", str (socket.gethostbyname(socket.gethostname())), wrapperFile, "-p", str (guess)])
    
    client_socket, address = server_socket.accept()
    
    client_socket.send (json.dumps ({"command": command, "iterations": iterations, "processes" : processes, "env": {"HOSTFILE": "hostfile_%04d" % guess }}).encode ())
    
    try:
        data = json.loads (client_socket.recv (512).decode ())
    except Exception as e:
        print (type (e), e)
    client_socket.close ()
    
    return data

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
            if not isinstance (arg, Argument):
                arg = Argument (*arg)
            self.variances.append (arg)
        self.uniques = kwargs.pop ("uniques", [])
        self.current = self.variances [:]
            
    def getSetupCommand (self, *args, **kwargs):
        command = [join (self.directory, self.setupCommandRoot)] + self.commandArgs
        for arg in self.variances + self.uniques:
            arg (command)
        return command
        
    def getCommand (self, *args, **kwargs):
        command = [join (self.directory, self.commandRoot)] + self.commandArgs
        arguments = [variance.copy (value = arg) for arg, variance in zip (args, self.variances [:len (args)])]
        arguments += self.variances [len (arguments):]
        for arg in arguments + self.uniques:
            arg (command, run = True)
        return command
        
    def getBaseTime (self):
        return timeCommand.delay (command = self.getCommand (*variances), setupCommand = self.getSetupCommand (*variances))
        
    basetime = property (getBaseTime)
        
    def calculateTimes (self, *args, **kwargs):
        times = {}
        if len (self.variances) == 0:
            return
            
        for variances in Argument.generate (*self.variances):
            if variances not in times:
                processes = 1
                threads = 1
                try:
                    Argument.setAll (self.variances, variances)
                except RuntimeError:
                    print ("Throwing out", variances)
                    continue
                    
                for arg in self.variances:
                    if arg.processes:
                        processes *= arg.value
                    if arg.threads:
                        threads *= arg.value

                for arg in self.uniques:
                    arg.setRandom ()
                times [variances] = timeCommand.delay (command = self.getCommand (), setupCommand = self.getSetupCommand (), processes = processes, threads = threads, commandRoot = self.commandRoot, **kwargs)

        return times

class Argument (object):
    """
    An object that can return a part of a command associated with varying a particular parameter
    """
    def __init__(self, command, extent = None, value = 1, **kwargs):
        super(Argument, self).__init__()
        self.command = command
        self.extent = extent
        self.lowerBound = kwargs.get ("lowerBound", 1)
        self.upperBound = kwargs.get ("upperBound", None)
        self.pastValues = []
        self.value = value
        self.kwargs = kwargs
        self.processes = kwargs.get ("processes", False)
        self.threads = kwargs.get ("threads", False)
        self.runOnly = kwargs.get ("runOnly", "")
        
    def getValue (self):
        return self._value
        
    def setValue (self, value):
        if (value < self.lowerBound) or (self.upperBound is not None and value > self.upperBound):
            raise RuntimeError ("Can't have a value below %s or above %s" % (str (self.lowerBound), str (self.upperBound)))
        if value not in self.pastValues:
            self.pastValues.append (value)
        self._value = value
        
    value = property (getValue, setValue)
        
    def setRandom (self):
        self.value += 1
        return self
        
    def __lt__ (self, other):
        if isinstance (other, Argument):
            return self.value - other.value < 0
        else: return self.value - other < 0
        
    def __gt__ (self, other):
        if isinstance (other, Argument):
            return self.value - other.value > 0
        else: return self.value - other > 0
        
    def __eq__ (self, other):
        if isinstance (other, Argument):
            return self.value == other.value
        else: return self.value == other
        
    def __mul__ (self, scalar):
        if isinstance (scalar, Argument):
            return CompositeArgument (self, scalar, operator.mul)
        x = copy.deepcopy (self)
        x.value = self.value * scalar
        return x
                
    def __add__ (self, scalar):
        if isinstance (scalar, Argument):
            return CompositeArgument (self, scalar, operator.add)
        x = copy.deepcopy (self)
        x.value = self.value + scalar
        return x
        
    def __div__ (self, scalar):
        if isinstance (scalar, Argument):
            return CompositeArgument (self, scalar, operator.div)
        x = copy.deepcopy (self)
        x.value = self.value / scalar
        return x
            
    def __sub__ (self, scalar):
        if isinstance (scalar, Argument):
            return CompositeArgument (self, scalar, operator.sub)
        x = copy.deepcopy (self)
        x.value = self.value - scalar
        return x
                
    def __call__ (self, fullCommand = None, run = False):
        if fullCommand is None:
            fullCommand = []
        if self.kwargs.get ("prepend", False):
            if run:
                for subcommand in (self.runOnly).split () [::-1]:
                    fullCommand.insert (0, subcommand)
            for subcommand in (self.command % self.value).split () [::-1]:
                fullCommand.insert (0, subcommand)
        else:
            for subcommand in (self.command % self.value).split ():
                fullCommand.append (subcommand)
            if run:
                for subcommand in (self.runOnly).split ():
                    fullCommand.append (subcommand)
        return fullCommand
        
    def __hash__ (self):
        return hash (self.command) + hash (self.value)
        
    def __str__ (self):
        return " ".join (self ())
    
    def __repr__ (self):
        return "<%s>" % str (self)
        
    @staticmethod
    def generate (*args):
        if len(args) > 1:
            for i in args [0].extent:
                for rest in Argument.generate (*args [1:]):
                    yield (i,) + rest
        else:
            if args [0].extent is None:
                yield (args [0].value,)
            else:
                for i in args [0].extent: yield (i,)
            
    @staticmethod
    def setAll (args, values):
        for arg, value in zip (args, values):
            if arg.extent is not None:
                arg.value = value

class CompositeArgument (Argument):
    def __init__ (self, argument1, argument2, operator = operator.add):
        self.__dict__ = copy.deepcopy (argument1.__dict__)
        self.argument1 = argument1
        self.argument2 = argument2
        self.operator = operator
        
    @property
    def value (self):
        return self.operator (self.argument1.value, self.argument2.value)
    