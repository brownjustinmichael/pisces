from datetime import datetime
from subprocess import call
import json

from sys import argv
import socket

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect(("localhost", int (argv [2])))

inputs = json.loads (client_socket.recv (512).decode ())

times = []

for i in range (inputs ["iterations"]):
    startTime = datetime.now ()
    print ("Calling", inputs ["command"])
    call (inputs ["command"])
    times.append ((datetime.now() - startTime).total_seconds ())

results = {}
results ["dt"] = sum (times) / len (times)

client_socket.send(json.dumps (results).encode ())
client_socket.close ()