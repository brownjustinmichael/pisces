from datetime import datetime
from subprocess import call
import json

from sys import argv
import socket

from optparse import OptionParser

def get_total_seconds(td): return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 1e6) / 1e6


parser = OptionParser()
parser.add_option("-p", "--port", dest="port")
parser.add_option("-a", "--address", dest = "address")
(options, args) = parser.parse_args()

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

print (options.address, options.port)
client_socket.connect((options.address, int (options.port)))

inputs = json.loads (client_socket.recv (512).decode ())

times = []

for i in range (inputs ["iterations"]):
    startTime = datetime.now ()
    print ("Calling", inputs ["command"])
    call (inputs ["command"])
    times.append (get_total_seconds ((datetime.now() - startTime)))

results = {}
results ["dt"] = sum (times) / len (times)

client_socket.send(json.dumps (results).encode ())
client_socket.close ()