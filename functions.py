import os
import sys
from pyDOE import lhs 
import numpy as np
import random
import pandas as pd

# make base directories for the analysis
def makeTotalDir(totaldir):
    try:
        os.makedirs(totaldir)
    except:
        pass

    try:
        os.makedirs(totaldir + "points/")
    except:
        pass

    try:
        os.makedirs(totaldir + "QVir_Analysis/")
    except:
        pass

# gets list of directories in the target directory
def getDirs(dir):
    startdir = os.getcwd()
    os.chdir(dir+"points")
    files = os.listdir(".")

    directories = []
    for item in files:
        if os.path.isdir(item) and item != 'QVir_Analysis': directories.append(item)

    os.chdir(startdir)

    
    directories.sort()  # to make sure the directories are sorted in case it doesnt complete
    return directories

# to see what ran last if runs dont complete
def update(cmd):
    outfile = open('last-ran.txt','a')
    outfile.write(cmd+'\n')
    outfile.close()

    print(cmd)

# runs the xml in jetscape
def runxml(xmlname):
    runjetscapecmd = "./runJetscape " + xmlname
    update(runjetscapecmd)
    os.system(runjetscapecmd)

# writing set of xmls to a file to be run later
def writeXmls(xmls,baseDir):
    outfile = open(baseDir+"/xml/cmds.txt","w")
    outfile.writelines(xml+"\n" for xml in xmls)
    outfile.close()

# reading xmls to run
def readXmls(baseDir):
    infile = open(baseDir+"/xml/cmds.txt","r")
    lines = infile.readlines()
    return lines

# sorting list of strings by int value
def intSort(dirs):
    dirs = [int(x) for x in dirs]
    dirs.sort()
    dirs = [str(x) for x in dirs]
    return dirs

# gets list of directories in the target directory
def concatDats(dir):
    startdir = os.getcwd()
    os.chdir(dir+"/dat")
    dats = os.listdir(".")

    cmd = "cat "
    for dat in dats:
        cmd = cmd + dat + " "
    cmd = cmd + " > " + dats[0].rstrip(".dat.gz") + "0.dat.gz"
    os.system(cmd)
    
    for dat in dats:
        os.remove(dat)

    os.chdir(startdir)