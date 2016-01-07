import os
import abc
import yaml
import subprocess
import collections

def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d

def Configuration(file_name=None, default=None, **kwargs):
    """
    Returns a dictionary setup object for use with code objects

    :param file: The file from which the configuration should be loaded
    :param kwargs: Any additional parameters that should be added can be added as kwargs
    :return: The configuration dictionary
    :rtype: `dict`
    """
    # Load the default configuration
    if default is not None:
        tmp = yaml.load(open(default))
    else:
        tmp = yaml.load(open(os.path.join(os.path.dirname(__file__), "../src/defaults.yaml")))

    # Load any additional keys from the given file
    if file_name is not None:
        tmp = update(tmp, yaml.load(open(file_name)))

    # Update any additional arguments provided to the function
    tmp.update (kwargs)

    # Add some defaults if they haven't been specified
    if "np" not in tmp:
        tmp ["np"] = 1
    if "wd" not in tmp:
        tmp ["wd"] = "."

    return tmp

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

    @property
    def np(self):
        return self.config ["np"]
    
    @property
    def wd(self):
        return self.config ["wd"]
    
    @abc.abstractmethod
    def setup(self):
        """
        Generates any files needed for the code to execute.
        """
        pass

    @abc.abstractmethod
    def call(self):
        """
        Returns a list of the executable arguments to run the code

        :return: A list of executable arguments that can be run by subprocess.call
        :rtype: `list`
        """
        pass

class ISCES(Code):
    """
    A class that translates from configuration to an executable form for ISCES
    """
    def __init__(self, config, config_file="config.yaml"):
        super(ISCES, self).__init__(config)
        self._config_file = config_file
        
    def setup(self):
        """
        Generates the config file containing the parameters for the run and the initial state of the system
        """
        cwd = os.getcwd()
        os.chdir(self.wd)
        with open (self._config_file, "w") as stream:
            stream.write (yaml.dump (self.config))
        # Check that the input and output directories exist and make them if they don't
        if not os.path.isdir(self.config ["root"] + self.config ["input"] ["directory"]):
            os.makedirs(self.config ["root"] + self.config ["input"] ["directory"])
        if not os.path.isdir(self.config ["root"] + self.config ["output"] ["directory"]):
            os.makedirs(self.config ["root"] + self.config ["output"] ["directory"])

        subprocess.call(["mpiexec", "-np", str(self.np), os.path.join(os.path.dirname(__file__), "../run/isces_init"), self._config_file])

        os.chdir(self.wd)

    def call(self):
        """
        Returns the ISCES executable.
        """
        return ["mpiexec", "-np", str(self.np), os.path.join(os.path.dirname(__file__), "../run/isces"), self._config_file]
        
class LauncherRegistry(type):
    registry = {}

    @staticmethod
    def register_class(target_class):
        LauncherRegistry.registry[target_class.__name__] = target_class

    def __new__(cls, clsname, bases, attrs):
        newclass = super(LauncherRegistry, cls).__new__(cls, clsname, bases, attrs)
        LauncherRegistry.register_class(newclass)  # here is your register function
        return newclass

class Launcher(object, metaclass=LauncherRegistry):
    """
    An "abstract" class that takes a Code object and executes it.

    If the base class is used, the code is executed directly without going through Celery or Torque, for example.
    """
    def __init__(self, code):
        super(Launcher, self).__init__()
        self.code = code

    def launch(self, init=True):
        if init:
            self.code.setup()
        subprocess.call(self.code.call ())
