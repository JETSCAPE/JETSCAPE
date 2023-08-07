#makes set of bayesian design points for pp collisions

import os
import multiprocessing as mp
from pickle import TRUE
from pyDOE import lhs 
from datetime import date
from functions import *
import pandas as pd
from operator import itemgetter
import time

# option reading
system = "LHC"
RHIC = False
Drun = False
Lrun = False
Erun = False
reading = False
rerunning = False
appending = False
softOnly = False
design = []
startdir = ""
for i, option in enumerate(sys.argv):
    if "RHIC" in option:
        system = "RHIC"
        RHIC = True
    if "Drun" in option:
        system = "Drun"
        Drun = True
        softOnly = True
    if "Lrun" in option:
        system = "Lrun"
        Lrun = True
        softOnly = True
    if "Erun" in option:
        system = "Erun"
        Erun = True
        softOnly = True
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

# date and file initiation
todays_date = date.today()
totaldir = "/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/runs/" + system + "-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"
if rerunning: 
    totaldir = startdir
    design = pd.read_csv(totaldir+"QVir_Analysis/parameters.txt")
    print("Using: "+totaldir+"QVir_Analysis/parameters.txt")

# function to make xmls
def makexml(bound, parameters, baseDir, xmltemplate, ECM):
    events  = 20000
    # extra events in first bin
    if RHIC:
        events = events*6
    if Drun:
        events = 20000
    if Lrun:
        events = 20000
    if Erun:
        events = 100000
    if bound[0] == 0:
        events = events*2

    # appending option and filename
    ogfile = baseDir + "/dat/PP_Bin" + str(bound[0])+ "_" + str(bound[1])
    if appending:
        filename = ogfile+"temp"
    else:
        filename = ogfile

    # building lines for xml file
    eventLine = "  <nEvents> " + str(events) + " </nEvents>"
    lowerLine = "      <pTHatMin>" + str(bound[0]) + "</pTHatMin>\n"
    upperLine = "      <pTHatMax>" + str(bound[1]) + "</pTHatMax>\n"
    fileLine = "  <outputFilename>" + filename + "</outputFilename>\n"
    
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
        if "<nEvents>" in line: 
            newlines.append(eventLine)
            continue
        elif "<pTHatMin>" in line:
            newlines.append(lowerLine)
            continue
        elif "<pTHatMax>" in line:
            newlines.append(upperLine)
            continue
        elif "<outputFilename>" in line:
            newlines.append(fileLine)
            continue

        writeold = True
        for paramline in paramlines:
            if paramline["name"] in line: 
                newlines.append(paramline["line"])
                writeold = False
                break

        if writeold: newlines.append(line)

    xmlname = baseDir+"/xml/pp_"+ECM+"_"+str(bound[0])+"_"+str(bound[1])+".xml"
    xml = open(xmlname,'w')
    xml.writelines(newlines)
    xml.close()

    output = {'xml': xmlname, 'dat': ogfile}
    return output

# making directory for the runs
def makeDir(index):
    # setting names
    baseDir = totaldir + str(index)
    plotsDir = baseDir + "/plots"
    rootDir = baseDir + "/root"
    datDir = baseDir + "/dat"
    xmlDir = baseDir + "/xml"

    #making dirs
    try:
        os.makedirs(plotsDir)
    except:
        pass
    
    try:
        os.makedirs(datDir)
    except:
        pass

    try:
        os.makedirs(rootDir)
    except:
        pass

    try:
        os.makedirs(xmlDir)
    except:
        pass

    #returning base directory name
    return baseDir

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Start of main method

# Creating total directory for use
try:
    os.makedirs(totaldir)
except:
    pass

try:
    os.makedirs(totaldir + "QVir_Analysis/")
except:
    pass

# pTHat bounds and xml based on RHIC vs LHC
intervals = []
ECM = ""
if RHIC:
    intervals = [0, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70]
    ECM = "200"
elif Drun or Erun:
    intervals = range(0,48+1) # need to add one more to the range than bins desired
    ECM = "5020"
elif Lrun:
    intervals = range(0,48+1) # need to add one more to the range than bins desired
    ECM = "13000"
else:
    intervals = [0, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 230, 250, 270, 290, 310, 330, 350, 400, 450, 500, 550, 600, 1000]
    ECM = "2760"

xmlname = "/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/config/jetscape_user_pp"+system+".xml"

xmltemplate = open(xmlname,"r")
xmllines = xmltemplate.readlines()

# making design points
nsamples = 48
if not reading: design = createPandaDesign(nsamples)

# pthat bounds formed into pairs
pTHatBounds = []
for i in range(len(intervals)-1):
    if softOnly:
        pTHatBounds.append((0,-i-1)) # setting all soft bins upper bound to negative
    elif i is 0 and intervals[i] == 0:
        pTHatBounds.append((intervals[i],-1))
    else:
        pTHatBounds.append((intervals[i],intervals[i+1]))

# Changing to build directory
os.chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build")

pool = mp.Pool(len(pTHatBounds))

# Loop over design points and make parameter file
design.to_csv(totaldir+'QVir_Analysis/parameters.txt')
print(design)

for i in range(len(design)):    
    # Making xmls
    baseDir = makeDir(i)
    output = pool.starmap(makexml, [(bound, design.loc[[i]], baseDir, xmllines, ECM) for bound in pTHatBounds]) # star map to each set of bounds

    # Running jetscape for them
    xmls, dats = zip(*map(itemgetter('xml', 'dat'), output)) 
    pool.map(runxml, xmls)

    # Handling appending runs
    if appending:
        for dat in dats:
            appendcommand = "cat " + dat + "temp.dat >> " + dat + ".dat"
            os.system(appendcommand)
            os.system('rm '+dat+'temp.dat')

    # concatonating all soft bins together
    if softOnly:
        cmd = "cat "
        for dat in dats:
            cmd = cmd + dat + ".dat.gz "
        cmd = cmd + " > " + dats[-1] + "0.dat.gz"
        os.system(cmd)
        
        for dat in dats:
            os.remove(dat+".dat.gz")

pool.close()