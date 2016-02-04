import os
import yaml
import abc
import subprocess
import shutil
import glob

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
        self.process = None

    def launch(self, *args, from_dump=True, **kwargs):
        if self.process is not None:
            raise RuntimeError("Already running")
        if from_dump:
            self.code.resume()
        else:
            self.code.setup()

        self._launch(*args, **kwargs)

    def _launch(self, *args, **kwargs):
        self.process = subprocess.Popen(self.code.call ())

    def wait(self, *args, record=True, session=None, **kwargs):
        if self.process is None:
            raise RuntimeError("Not yet running")
        self._wait(*args, **kwargs)
        if record:
            self.code.record(session)
        self.process = None

    def _wait(self, *args, **kwargs):
        self.process.wait()

