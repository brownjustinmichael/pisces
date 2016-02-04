from celery import Celery
import json
import operator

from . import celeryconfig
from .code import CodeRegistry
from .launch import Launcher, LauncherRegistry

app = Celery('pisces', backend = 'amqp')
app.config_from_object (celeryconfig)

@app.task
def launch(class_name, config, launcher="Launcher"):
    code = CodeRegistry.registry [class_name] (config)
    launcher = LauncherRegistry.registry [launcher] (code)
    return launcher._launch()

class CeleryLauncher(Launcher):
    """
    A class designed to launch code execution through celery
    """
    def __init__(self, code, sub_launcher="Launcher"):
        super().__init__(code)
        self.code = code
        self.sub_launcher = sub_launcher

    def _launch(self, *args, **kwargs):
        self.process = launch.delay(self.code.__class__.__name__, self.code.config, launcher=self.sub_launcher)

    def _wait(self, *args, **kwargs):
        self.process.get ()
