from datetime import datetime
from subprocess import call
import json

from sys import argv
import socket

def get_total_seconds(td): return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 1e6) / 1e6

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect(("localhost", int (argv [2])))

inputs = json.loads (client_socket.recv (512).decode ())

times = []

for i in range (inputs ["iterations"]):
    startTime = datetime.now ()
    print ("Calling", inputs ["command"])
    call (inputs ["command"])
    times.append (get_total_second ((datetime.now() - startTime)))

results = {}
results ["dt"] = sum (times) / len (times)

client_socket.send(json.dumps (results).encode ())
client_socket.close ()