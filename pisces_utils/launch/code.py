import os
import yaml
import abc
import subprocess
import shutil
import glob
import re
import datetime

try:
    from .. import db
except ImportError as e:
    db = None
    print("Could not import database (%s)" % str(e))

class CodeRegistry(type):
    registry = {}

    @staticmethod
    def register_class(target_class):
        CodeRegistry.registry[target_class.__name__] = target_class

    def __new__(cls, clsname, bases, attrs):
        newclass = super(CodeRegistry, cls).__new__(cls, clsname, bases, attrs)
        CodeRegistry.register_class(newclass)  # here is your register function
        return newclass

class Code(object, metaclass=CodeRegistry):
    """
    An abstract class describing the way to translate from a configuration to code execution
    """
    def __init__(self, config):
        super(Code, self).__init__()
        self.config = config

    @classmethod
    def from_simulation_entry(cls, simulation_entry):
        return cls(simulation_entry.get_config())

    @property
    def np(self):
        return self.config["np"]

    @np.setter
    def np(self, value):
        self.config["np"] = value
    
    @property
    def wd(self):
        return self.config ["wd"]

    @property
    def date(self):
        """Return the last modified time of the executable"""
        return datetime.datetime.fromtimestamp (os.path.getmtime(self.executable()))
    
    @abc.abstractmethod
    def setup(self):
        """
        Generates any files needed for the code to execute.
        """
        pass

    @abc.abstractmethod
    def resume(self):
        """
        Generates any files needed for the code to resume a previous run.
        """
        pass

    def call(self, x=None, mpi="mpiexec", debug=2, **kwargs):
        """
        Returns the ISCES executable.
        """
        if x is None:
            x = []
        call = [mpi]
        for arg in x:
            call += ["-x", arg]
        for arg in kwargs:
            call += ["-" + arg, str(kwargs[arg])]
        return call + ["-n", str(self.np), self.executable(), self._config_file, "-D%i" % debug]

    @abc.abstractmethod
    def executable(self):
        """
        Returns a list of the executable arguments to run the code

        :return: A list of executable arguments that can be run by subprocess.call
        :rtype: `list`
        """
        pass

    @abc.abstractmethod
    def record(self):
        """
        Stores the run in a database
        """
        pass

class PISCES(Code):
    """
    A class that translates from configuration to an executable form for ISCES
    """
    def __init__(self, config, config_file="config.yaml", init=None):
        super(PISCES, self).__init__(config)
        self._config_file = config_file
        if init is not None:
            config["init"] = init
        else:
            if "init" not in config:
                config["init"] = "init"

    @property
    def threads(self):
        return self.config ["parallel"] ["maxthreads"]
        
    def setup(self, init=True):
        """
        Generates the config file containing the parameters for the run and the initial state of the system
        """
        cwd = os.getcwd()
        print("Setting up in ", self.wd)
        if not os.path.isdir(self.wd):
            os.makedirs(self.wd)
        os.chdir(self.wd)
        with open (self._config_file, "w") as stream:
            stream.write (yaml.dump (self.config))
        # Check that the input and output directories exist and make them if they don't
        for dirname in [self.config ["input"] ["directory"], self.config ["output"] ["directory"], self.config ["dump"] ["directory"]]:
            if not os.path.isdir(self.config ["root"] + dirname):
                os.makedirs(self.config ["root"] + dirname)

        if init:
            subprocess.call(["mpiexec", "-np", str(self.np), os.path.join(os.path.dirname(__file__), "../../run/", self.config ["init"]), self._config_file])

        os.chdir(self.wd)

    def resume(self):
        cwd = os.getcwd()
        os.chdir(self.wd)

        try:
            os.makedirs(os.path.join(self.config ["dump"] ["directory"], "previous_%02i" % self.config ["output"] ["number"]))
        except FileExistsError:
            pass
        for file in glob.glob(os.path.join(self.config ["dump"] ["directory"], self.config ["dump"] ["file"].replace("%02i", "*")) + ".cdf"):
            shutil.copy(file, os.path.join(self.config ["dump"] ["directory"], "previous_%02i" % self.config ["output"] ["number"]))
        self.config ["input"] ["file"] = self.config ["dump"] ["file"]
        self.config ["input"] ["directory"] = self.config ["dump"] ["directory"]
        self.config ["output"] ["number"] += 1
        self.setup(init=False)

        os.chdir(self.wd)


    def executable(self):
        """
        Returns the ISCES executable.
        """
        return os.path.join(os.path.dirname(__file__), "../../run/pisces")

    def record(self, session=None):
        if db is None:
            raise RuntimeError("Could not record as the database backend is not functional")
        if session is None:
            session = db.Session()

        filename = os.path.join(self.config["root"], self.config["output"]["directory"], self.config["output"]["stat"]["file"])
        filename = re.sub("\%*\d*[id]", "[0-9]*", filename)
        entry = None

        entry = db.SimulationEntry.from_config(session, **self.config)
        entry.date = self.date
        session.add (entry)

        for filename in glob.glob(filename + ".*"):
            steps = db.StepEntry.from_file (session, filename, sim = entry)
            for step in steps:
                session.add (step)

        session.commit ()
        return entry
        