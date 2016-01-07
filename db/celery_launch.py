from celery import Celery
import celeryconfig
import json
import operator

from launch import Launcher, Registry

app = Celery('pisces_timer', backend = 'amqp')
app.config_from_object (celeryconfig)

@app.task
def launch(self, class_name, config):
    code = Registry.registry [class_name] (config)
    launcher = Launcher (code)
    return launcher.launch()

class CeleryLauncher(Launcher):
    """
    A class designed to launch code execution through celery
    """
    def __init__(self, code):
        super(Launcher, self).__init__()
        self.code = code

    @app.task
    def launch(self, init=True):
        launch.delay(self.code.__name__, self.code.config)