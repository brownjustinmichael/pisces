import sys
import argparse

parser = argparse.ArgumentParser ()
parser.add_argument ('pbs_file')
parser.add_argument ('--ppn', default = None, type = int)

namespace = parser.parse_args ()

file = open (sys.argv [1], "r")

hosts = []

if namespace.ppn is None:
    for line in file:
        if line [:-1] not in hosts and line != "" and line != "\n":
            hosts.append (line [:-1])
else:
    for i, line in enumerate (file):
        if i % namespace.ppn == 0:
            hosts.append (line [:-1])

for host in hosts:
    print ("%s:1" % host)
