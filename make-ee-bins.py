#makes set of bayesian design points for ee collisions

import os
import multiprocessing as mp
import pandas as pd
from datetime import date
from functions import *

#option reading
reading = False
for i, option in enumerate(sys.argv):
    if "-d" in option:
        reading = True
        design = pd.read_csv(sys.argv[i+1])
    if "-r" in option:
        rerunning = True
        reading = True
        startdir = sys.argv[i+1]
    if "-a" in option:
        appending = True
        rerunning = True
        reading = True
        startdir = sys.argv[i+1]

#date and file initiation
todays_date = date.today()
totaldir = "/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/runs/LEP-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"

makeTotalDir(totaldir)

#methods for runnning
def binrun(index, parameters, xmltemplate):
    ##making directory for this designpoint
    baseDir = makeDir(index)

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

    ##running for the parameters set in the xml file
    runxml(xmlname)

    analysiscmd = "./ee-analysis-spectra " + baseDir
    os.system(analysiscmd)

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

#reading xml template
xmltemplate = open("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/config/jetscape_user.xml","r")
xmllines = xmltemplate.readlines()

##Changing to build directory
os.chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build")

# Loop over design points and make parameter file
design.to_csv(totaldir+'QVir_Analysis/parameters.txt')
print(design)

##Running jetscape for each set of parameters with multiprocessing
pool = mp.Pool(48)
pool.starmap(binrun, [(i,design.loc[[i]], xmllines) for i in range(len(design))])
pool.close()

##running total analysis
#os.system("./ee-comparison " + totaldir)