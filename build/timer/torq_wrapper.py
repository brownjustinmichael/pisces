from datetime import datetime
from subprocess import call
import json

from sys import argv
import socket

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect(("localhost", int (argv [2])))

inputs = json.loads (client_socket.recv (512).decode ())

print ("Writing to batch file")

batch_file = open ("batch_%04d.pbs" % int (argv [2]), "w")

batch_file.write ("#PBS -S /bin/bash\n")
batch_file.write ("#PBS -q normal\n")
batch_file.write ("#PBS -N pisces\n")

if (inputs ["processors"] > 16):
    ppn = 16
else:
    ppn = inputs ["processors"]
    
batch_file.write ("#PBS -l nodes=%d:ppn=%d\n" % (inputs ["processors"] / 16 + 1, ppn))
batch_file.write ("#PBS -l walltime=00:30:00\n")
batch_file.write ("cd $PBS_O_WORKDIR\n")
batch_file.write ("cp $PBS_NODEFILE .\n")

batch_file.write (" ".join (inputs ["command"]))
batch_file.write ("\n")

batch_file.close ()

times = []

for i in range (inputs ["iterations"]):
    startTime = datetime.now ()
    # call (["qsub", "batch_&04d.pbs" % int (argv [2])])
    times.append ((datetime.now() - startTime).total_seconds ())

results = {}
results ["dt"] = sum (times) / len (times)

client_socket.send(json.dumps (results).encode ())
client_socket.close ()