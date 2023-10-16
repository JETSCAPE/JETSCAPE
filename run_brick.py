#makes set of bayesian design points for pp collisions

from fileinput import filename
import os
import multiprocessing as mp
from pickle import TRUE
from datetime import date
from functions import *
import pandas as pd
from operator import itemgetter
import time
import xml.etree.ElementTree as ET

# option reading
softOnly = True
design = pd.read_csv("designs/default.txt")
startdir = ""

# date and file initiation
todays_date = date.today()
totaldir = "/scratch/user/cameron.parker/newJETSCAPE/gamma/JETSCAPE/runs/Brick-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"

# xml shenanigans
xmltemplate = "/scratch/user/cameron.parker/newJETSCAPE/gamma/JETSCAPE/config/jetscape_user_brick_hybrid_hadronization.xml"

# function to make xmls
def makexml(index, parameters, baseDir):
    # appending option and filename
    filename = baseDir + "/dat/Brick_Bin" + str(index)

    # xml declarations
    newxml = ET.parse(xmltemplate)
    root = newxml.getroot()

    # building lines for xml file
    for temp in root.iter('outputFilename'):
        temp.text = str(filename)
    
    for key in parameters.columns:
        for temp in root.iter(key):
            temp.text = str(parameters.iloc[-1][key])

    xmlname = baseDir+"/xml/brick_"+str(index)+".xml"
    newxml.write(xmlname)

    output = {'xml': xmlname, 'dat': filename}
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

# pthat bounds formed into pairs
indices = list(range(1))

# Changing to build directory
os.chdir("/scratch/user/cameron.parker/newJETSCAPE/gamma/JETSCAPE/build")

pool = mp.Pool(len(indices))

# Loop over design points and make parameter file
design.to_csv(totaldir+'QVir_Analysis/parameters.txt',index=False)
print(design)

for i in range(len(design)):    
    # Making xmls
    baseDir = makeDir(i)
    output = pool.starmap(makexml, [(index, design.loc[[i]], baseDir) for index in indices]) # star map to each set of bounds

    # Running jetscape for them
    xmls, dats = zip(*map(itemgetter('xml', 'dat'), output)) 
    writeXmls(xmls,baseDir)
    #continue
    pool.map(runxml, xmls)

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