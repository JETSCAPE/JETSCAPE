#makes set of bayesian design points for ee collisions

import os
import multiprocessing as mp
import pandas as pd
from datetime import date
from functions import *
import htcondor

#option reading
reading = False
blank = False
parallel = False
todays_date = date.today()
name = "LEP-" + str(todays_date.month) + "-" + str(todays_date.day)
for i, option in enumerate(sys.argv):
    if "-d" == option:
        reading = True
        design = pd.read_csv(sys.argv[i+1])
    if "-r" == option:
        rerunning = True
        reading = True
        startdir = sys.argv[i+1]
    if "-a" == option:
        appending = True
        rerunning = True
        reading = True
        startdir = sys.argv[i+1]
    if "-n" == option:
        name = sys.argv[i+1]
    if "--blank" == option:
        blank = True
    if "-p" == option:
        parallel = True

#date and file initiation
totaldir = "/data/rjfgroup/rjf01/cameron.parker/runs/" + name + "/"
makeTotalDir(totaldir)

# condor stuff
thisjob = xmljob
schedd = htcondor.Schedd()                   # get the Python representation of the scheduler

#methods for runnning
def binrun(index, parameters, xmltemplate):
    ##making directory for this designpoint
    baseDir = makeDir(index)

    # link QS to QO
    if 'QSfactor' in parameters.columns:
        parameters['QS'] = (2*parameters.lambdaQCD+0.05) + (parameters.Q0-(2*parameters.lambdaQCD+0.05))*parameters.QSfactor

    ##building lines for xml file
    fileLine = "    <outputFilename>" + baseDir + "/run" + "</outputFilename>\n"
    paramlines = []
    for key in parameters.columns:
        paramline = {}
        paramline["name"] = key
        if ":" in key:
            paramline["line"] = "\t\t\t"+key+" = "+str(parameters.iloc[-1][key])+"\n"
        else:
            paramline["line"] = "\t\t\t<"+key+">"+str(parameters.iloc[-1][key])+"</"+key+">\n"
        paramlines.append(paramline)

    newlines = []
    for line in xmltemplate:
        if "<outputFilename>" in line:
            newlines.append(fileLine)
            continue

        writeold = True
        for paramline in paramlines:
            if paramline["name"] in line: 
                newlines.append(paramline["line"])
                writeold = False
                break

        if writeold: newlines.append(line)

    xmlname = baseDir + "/config.xml"
    xml = open(xmlname,'w')
    xml.writelines(newlines)
    xml.close()

    if blank:
        return

    ##running for the parameters set in the xml file
    if parallel:
        xmlinput = [{"xml": xmlname}]
        thisjob["batch_name"] = name + "-" + str(index)
        submit_result = schedd.submit(thisjob, itemdata = iter(xmlinput))
        print("Submitted job for:",index)
    else:
        runxml(xmlname)

##making directory for the runs
def makeDir(index):
    ##setting names
    baseDir = totaldir + "points/" + str(index)

    ##making dirs
    try:
        os.makedirs(baseDir)
        os.makedirs(baseDir+"/plots")
    except:
        pass

    ##returning base directory name
    return baseDir

######################################################################################################################################
##start of main method

##parameters and bounds
nsamples = 500
if not reading: design = createPandaDesign(nsamples)
try: design.drop(['vir_factor','MultipartonInteractions:ecmPow','MultipartonInteractions:pT0Ref'], axis=1, inplace=True)
except: temp = 1
design = design.rename(columns={'ee_vir_factor': 'vir_factor'})
design.to_csv(totaldir+'QVir_Analysis/parameters.txt',index=False)
print(design)

#reading xml template
xmltemplate = open("/data/rjfgroup/rjf01/cameron.parker/builds/JETSCAPE/config/jetscape_user.xml","r")
xmllines = xmltemplate.readlines()

##Changing to build directory
os.chdir("/data/rjfgroup/rjf01/cameron.parker/builds/JETSCAPE/build")

##Running jetscape for each set of parameters with multiprocessing
for i in range(len(design)):
    binrun(i,design.loc[[i]], xmllines)