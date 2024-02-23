#makes set of bayesian design points for pp collisions

import os
import multiprocessing as mp
from pickle import TRUE
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
LHC13000 = False
reading = False
rerunning = False
appending = False
softOnly = False
gzip = True
ECM = ""
design = []
startdir = ""
todays_date = date.today()
name = system + "-" + str(todays_date.month) + "-" + str(todays_date.day)
rundata = True

for i, option in enumerate(sys.argv):
    if "RHIC" in option:
        system = "RHIC"
        RHIC = True
        ECM = "200"
    if "LHC13000" in option:
        system = "LHC13000"
        LHC13000 = True
        ECM = "LHC13000"
    if "Drun" in option:
        system = "Drun"
        Drun = True
        gzip = False
        softOnly = True
        ECM = "5020"
    if "Lrun" in option:
        system = "Lrun"
        Lrun = True
        gzip = False
        softOnly = True
        ECM = "7000"
    if "Erun" in option:
        system = "Erun"
        Erun = True
        softOnly = True
        ECM = "5020"
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
    if "-n" in option:
        name = sys.argv[i+1]
    if "--blank" in option:
        rundata = False
        

# date and file initiation
totaldir = "/scratch/user/cameron.parker/projects/JETSCAPE/runs/" + name + "/"
if rerunning: 
    totaldir = startdir
    design = pd.read_csv(totaldir+"QVir_Analysis/parameters.txt")
    print("Using: "+totaldir+"QVir_Analysis/parameters.txt")

# function to make xmls
def makexml(bound, parameters, baseDir, xmltemplate, ECM):
    # appending option and filename
    ogfile = baseDir + "/dat/PP_Bin" + str(bound[0])+ "_" + str(bound[1])
    if appending:
        filename = ogfile+"temp"
    else:
        filename = ogfile

    # building lines for xml file
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
        if "<pTHatMin>" in line:
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
    baseDir = totaldir + "points/" + str(index)
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
makeTotalDir(totaldir)

# pTHat bounds and xml based system being run
intervals = []
if RHIC:
    intervals = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70]
else:
    intervals = [0, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 400, 450, 500, 550, 600, 1000]

if softOnly:
    intervals = range(0,48+1) # need to add one more to the range than bins desired

xmlname = "/scratch/user/cameron.parker/projects/JETSCAPE/config/jetscape_user_pp"+system+".xml"

xmltemplate = open(xmlname,"r")
xmllines = xmltemplate.readlines()

# making design points
nsamples = 500
if not reading: design = createPandaDesign(nsamples)

# pthat bounds formed into pairs
pTHatBounds = []
for i in range(len(intervals)-1):
    if softOnly:
        pTHatBounds.append((0,-i-1)) # setting all soft bins upper bound to negative
    elif intervals[i] == 0:
        pTHatBounds.append((intervals[i],-i-1))
    else:
        pTHatBounds.append((intervals[i],intervals[i+1]))

print(len(pTHatBounds))

# Changing to build directory
os.chdir("/scratch/user/cameron.parker/projects/JETSCAPE/build")

pool = mp.Pool(len(pTHatBounds))

# Loop over design points and make parameter file
design = design.loc[:, ~design.columns.str.contains('^Unnamed')]
try: design.drop(columns=['ee_vir_factor'],inplace=True)
except: temp = 1
design.to_csv(totaldir+'QVir_Analysis/parameters.txt',index=False)
print(design)

for i in range(len(design)):    
    # Making xmls
    baseDir = makeDir(i)
    output = pool.starmap(makexml, [(bound, design.loc[[i]], baseDir, xmllines, ECM) for bound in pTHatBounds]) # star map to each set of bounds
    
    # Running jetscape for them
    xmls, dats = zip(*map(itemgetter('xml', 'dat'), output)) 
    writeXmls(xmls,baseDir)
    if rundata: pool.map(runxml, xmls)
    zipxmls(baseDir)

    if not rundata: continue
    # Handling appending runs
    if appending:
        for dat in dats:
            appendcommand = "cat " + dat + "temp.dat >> " + dat + ".dat"
            os.system(appendcommand)
            os.system('rm '+dat+'temp.dat')

    # concatonating all soft bins together
    if "Bin0" in dats[1]:
        softCombine(baseDir)

pool.close()