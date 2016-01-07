from celery import Celery
import db.celeryconfig as celeryconfig
import json
import operator

from db.launch import Launcher, CodeRegistry, LauncherRegistry

app = Celery('pisces', backend = 'amqp')
app.config_from_object (celeryconfig)

@app.task
def launch(self, class_name, config, launcher="Launcher"):
    code = CodeRegistry.registry [class_name] (config)
    launcher = LauncherRegistry.registry [launcher] (code)
    return launcher.launch()

class CeleryLauncher(Launcher):
    """
    A class designed to launch code execution through celery
    """
    def __init__(self, code, sub_launcher="Launcher"):
        super(Launcher, self).__init__()
        self.code = code
        self.sub_launcher = sub_launcher

    @app.task
    def launch(self, init=True):
        launch.delay(self.code.__name__, self.code.config, launcher=self.sub_launcher)