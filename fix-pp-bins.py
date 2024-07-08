from functions import *
import os
import htcondor
import re
import sys

# option reading
runAll = False
for i, option in enumerate(sys.argv):
    if "-all" == option:
        runAll = True

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
directories = intSort(directories)

# condor initializations
testjob = xmljob
schedd = htcondor.Schedd()

# directory loop
for directory in directories:
    thisDir = analysisDir + "points/" + directory
    xmlsToRun = []

    xmls = readXmls(thisDir)
    xmls.sort()

    dats = os.listdir(thisDir+"/dat")
    dats.sort()

    # finding out what needs to be rerun
    for xml in xmls:
        if runAll:
            xmlsToRun.append({"xml": xml})
        else:
            run = True
            startBound = xml.split("_")[-2]

            for dat in dats:
                if("nfs" in dat): continue

                datBound = dat.split("_")[1].split("n")[1]
                if startBound == datBound:
                    run = False
                    break
            
            if run:
                xmlsToRun.append({"xml": xml})
            
    # submission
    if len(xmlsToRun) > 0:
        testjob["batch_name"] = analysisDir.split('/')[-2] + "-" + directory
        submit_result = schedd.submit(testjob, itemdata = iter(xmlsToRun))
        print(len(xmlsToRun),"fixes sumbitted for", directory)
    #print(xmlsToRun)