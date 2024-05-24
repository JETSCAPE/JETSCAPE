from functions import *
import os
import htcondor
import re

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
directories = intSort(directories)

# condor initializations
testjob = xmljob
schedd = htcondor.Schedd()

# directory loop
for directory in directories:
    run = True
    thisDir = analysisDir + "points/" + directory

    files = os.listdir(thisDir)
    files.sort()
    for file in files:
        if "run.dat.gz" in file:
            run = False

    if run:
        xmlname = [{"xml": thisDir+"/config.xml"}]
        testjob["batch_name"] = analysisDir.split('/')[-2] + "-" + directory
        submit_result = schedd.submit(testjob, itemdata = iter(xmlname))
        print("Fix sumbitted for", directory)