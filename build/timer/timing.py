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
def timeCommand (command, setupCommand = None, iterations = 1, wrapperFile = "wrapper.py", processors = 1):
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

    Popen (["python", wrapperFile, "-p", str (guess)])
    
    client_socket, address = server_socket.accept()
    
    client_socket.send (json.dumps ({"command": command, "iterations": iterations, "processors" : processors}).encode ())
    
    data = json.loads (client_socket.recv (512).decode ())
    client_socket.close ()
    
    return float (data ["dt"])

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
        for arg in self.variances + self.uniques:
            arg (command)
        return command
        
    def getBaseTime (self):
        return timeCommand.delay (command = self.getCommand (*variances), setupCommand = self.getSetupCommand (*variances))
        
    basetime = property (getBaseTime)
        
    def calculateTimes (self, *args, **kwargs):
        times = {}
        if len (self.variances) == 0:
            return
            
        print (self.variances)
        for variances in Argument.generate (*self.variances):
            if variances not in times:
                processors = 1
                for arg in self.variances:
                    if arg.threaded:
                        processors *= arg.value
                try:
                    for arg in self.uniques:
                        arg.setRandom ()
                    Argument.setAll (self.variances, variances)
                    times [variances] = timeCommand (command = self.getCommand (), setupCommand = self.getSetupCommand (), processors = processors, **kwargs)
                except RuntimeError:
                    print ("Throwing out", variances)
                    pass
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
        self.threaded = kwargs.get ("threaded", False)
        
    def getValue (self):
        return self._value
        
    def setValue (self, value):
        if (value < self.lowerBound) or (self.upperBound is not None and value > self.upperBound):
            raise RuntimeError ("Can't have a value below %s or above %s" % (str (self.lowerBound), str (self.upperBound)))
        if value not in self.pastValues:
            self.pastValues.append (value)
        self._value = value
        
    value = property (getValue, setValue)
        
    def copy (self, **kwargs):
        newkwargs = copy.copy (self.kwargs)
        newkwargs ["value"] = self.value
        newkwargs ["extent"] = self.extent
        for arg in kwargs:
            newkwargs [arg] = kwargs [arg]
        newObject = Argument (self.command, **newkwargs)
        return newObject
        
    def setRandom (self):
        self.value += 1
        return self
        
    def __lt__ (self, other):
        return self - other < 0
        
    def __gt__ (self, other):
        return self - other > 0
        
    def __eq__ (self, other):
        return self - other == 0
        
    def __mul__ (self, scalar):
        if isinstance (scalar, Argument):
            scalar = scalar.value
        return self.value * scalar
                
    def __add__ (self, scalar):
        if isinstance (scalar, Argument):
            scalar = scalar.value
        return self.value + scalar
        
    def __div__ (self, scalar):
        if isinstance (scalar, Argument):
            scalar = scalar.value
        return self.value / scalar
            
    def __sub__ (self, scalar):
        if isinstance (scalar, Argument):
            scalar = scalar.value
        return self.value - scalar
                
    def __call__ (self, fullCommand = None):
        if fullCommand is None:
            fullCommand = []
        if self.kwargs.get ("prepend", False):
            for subcommand in (self.command % self.value).split () [::-1]:
                fullCommand.insert (0, subcommand)
        else:
            for subcommand in (self.command % self.value).split ():
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
            for i in args [0].extent: yield (i,)
            
    @staticmethod
    def setAll (args, values):
        for arg, value in zip (args, values):
            arg.value = value

class CompositeArgument (object):
    def __init__ (self, argument1, argument2, operator = operator.add):
        self.argument1 = argument1
        self.argument2 = argument2
        self.operator = operator
        
    def getValue (self):
        return self.operator (self.argument1, self.argument2)
        
    value = property (getValue)
        
    def __add__ (self, other):
        return self.value + other.value
        
    def __sub__ (self, other):
        return self.value - other.value
        
    def __mul__ (self, other):
        return self.value * other.value
        
    def __div__ (self, other):
        return self.value / other.value
        
    def __gt__ (self, other):
        return self - other > 0
        
    def __lt__ (self, other):
        return self - other < 0
        
    def __eq__ (self, other):
        return self - other == 0
    
    