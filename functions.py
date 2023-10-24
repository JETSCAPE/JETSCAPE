from math import sqrt
import os
import sys
from pyDOE import lhs 
import numpy as np
import random
import pandas as pd
import math

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
        os.makedirs(totaldir + "analysis/")
    except:
        pass

# gets list of directories in the target directory
def getDirs(dir):
    startdir = os.getcwd()
    os.chdir(dir+"points")
    files = os.listdir(".")

    directories = []
    for item in files:
        if os.path.isdir(item) and item != 'analysis': directories.append(item)

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

class parton:
    def __init__(self, index, PID, status, E, Px, Py, Pz):
        self.index = index
        self.PID = PID
        self.status = status
        self.E = E
        self.Px = Px
        self.Py = Py
        self.Pz = Pz
        self.pT = sqrt(Px**2 + Py**2)

    def __init__(self, line):
        params = line.split(' ')
        if len(params) < 7:
            self.index = -1
            self.PID = -1
            self.status = -1
            self.E = -1
            self.Px = -1
            self.Py = -1
            self.Pz = -1
            self.pT = -1
            return

        self.index = int(params[0])
        self.PID = int(params[1])
        self.status = int(params[2])
        self.E = float(params[3])
        self.Px = float(params[4])
        self.Py = float(params[5])
        self.Pz = float(params[6])
        self.pT = sqrt(self.Px**2 + self.Py**2)

    def __str__(self):
        return f"{self.PID} {self.E} {self.pT}"
    
##chance for photon to be absorbed from https:##arxiv.org/abs/hep-ph/9405309
def absFactor1(pVec, T):
    alpha = 1./137.
    alphaS = 0.3
    p = pVec.P()
    return 2.0*(5.*math.pi/9.)*(alpha*alphaS*T*T/p)*math.log(0.2317*p/(alphaS*T))

def absFactor2(pVec, T):
    ## Constants
    alpha = 1.0 / 137.0 
    hbc = 0.1973 
    
    ## Calculate alpha_s at temperature 'temp'
    alphsT = 6.0 * math.pi / (27.0 * math.log(T / 0.022)) 
    gsT = math.sqrt(alphsT * 4.0 * math.pi) 
    
    ## Calculate x and exponential term
    p = pVec.P() 
    x = pVec.P() / T 
    expo = math.exp(-x) 
    fermi = expo / (1.0 + expo) 
    
    ## Calculate prfph
    prfph = alpha * alphsT * T**2 * (5.0 / 9.0) / (p * math.pi**2) 
    
    ## Calculate C22 and Cab
    C22 = 0.041 / x - 0.3615 + 1.01 * math.exp(-1.35 * x) 
    Cab = math.sqrt(1.5) * (0.548 / pow(x, 1.5) * math.log(12.28 + 1.0 / x) + 0.133 * x / math.sqrt(1.0 + x / 16.27)) 
    
    ## Calculate Ctot
    Ctot = 0.5 * math.log(2.0 * x) + C22 + Cab 
    
    ## Calculate dRd3p
    dRd3p = prfph * fermi * (math.log(math.sqrt(3.0) / gsT) + Ctot)
    return dRd3p * 4 * pow(math.pi,3) / expo