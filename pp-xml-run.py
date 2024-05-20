# running over points for pp collisions
# call with python3 pp-xml-run.py [DIRECTORY] [START] [END]

import os
from functions import *
import htcondor
import time
import string

totdir = sys.argv[1]
dirs = getDirs(totdir)

# job initialization
testjob = xmljob
schedd = htcondor.Schedd()                   # get the Python representation of the scheduler

# selecting range to run over
start = 0
finish = len(dirs)
dirs = intSort(dirs)

start = int(sys.argv[2])
if len(sys.argv) == 4:
    finish = int(sys.argv[3]) + 1

print("From " + str(start) + " to " + str(finish) + " in " + totdir + ":")

# looping over points
for dir in dirs[start:finish]:
    pointdir = totdir+"points/"+dir
    xmls = readXmls(pointdir)
    xmlinput = [{"xml": xml.rstrip()} for xml in xmls]
    testjob["batch_name"] = totdir.split('/')[-2] + "-" + dir
    time.sleep(3)
    submit_result = schedd.submit(testjob, itemdata = iter(xmlinput))
    print("Jobs sumbitted for:", dir)

    # concatonating all soft bins together
    #if "_0_" in xmls[0] and "_0_" in xmls[1]:
    #    dats = softCombine(pointdir)    