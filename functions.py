import os
import sys
from pyDOE import lhs 
import numpy as np
import random
import pandas as pd
from math import sqrt

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

# reading design file from an existing run
def readDesign(designfile):
    design = []

    # opening file
    file = open(designfile, 'r')
    lines = file.readlines()

    index = 0
    for line in lines:
        # kipping comments
        if '#' in line:
            continue
        
        # splitting line and converting to floats
        params = line.split()
        newparams = [index]
        for param in params:
            newparams.append(float(param))
        design.append(newparams)
        index = index + 1

    return design

# ranges of form [(start,end),(start,end)]
def createDesign(points, designparams):
    # initialization
    rawdesign = lhs(len(designparams), samples = points)
    index = 0
    design = []

    # Rescaling LHS parameters to appropriate ranges
    for point in rawdesign:
        newpoint = []
        newpoint.append({"index": index})
        for i, param in enumerate(designparams):
            newpoint.append({param["name"]: (point[i]*(param["range"][1]-param["range"][0]) + param["range"][0])})

        design.append(newpoint)
        index = index + 1

    return design

# ranges of form [(start,end),(start,end)]
def createPandaDesign(points):
    # design parameters
    designparams = [
        {"name": "Q0", "range": (0.95,3.0)},
        {"name": "vir_factor", "range": (0.1,1.0)},
        {"name": "lambdaQCD", "range": (0.1,0.4)},
        {"name": "part_prop", "range": (0.0,1.5)},
        {"name": "StringFlav:probStoUD", "range": (0.2,0.5)},
        {"name": "StringFlav:probQQtoQ", "range": (0.07,0.2)},
        {"name": "MultipartonInteractions:ecmPow", "range": (0.0,0.25)},
        {"name": "MultipartonInteractions:pT0Ref", "range": (0.5,2.5)},
        {"name": "pionWidthScale", "range": (0.5,2.0)},
        {"name": "kaonWidthScale", "range": (0.5,2.0)},
        {"name": "protonWidthScale", "range": (0.5,2.0)},
    ]

    # initialization
    rawdesign = lhs(len(designparams), samples = points)
    design = []

    # Rescaling LHS parameters to appropriate ranges
    for point in rawdesign:
        newpoint = []
        for i, param in enumerate(designparams):
            newpoint.append((point[i]*(param["range"][1]-param["range"][0]) + param["range"][0]))

        design.append(newpoint)

    df = pd.DataFrame(design,columns=[designparam["name"] for designparam in designparams])
    return df

# augment latin hypercube
def augmentDesign(points, bounds, n):
    m = len(points)
    nparams = len(points[0])-1
    temppoints = points

    # setting up intervals for hypercube cells
    intervals = []
    for bound in bounds:    # loop over parameters
        thisinterval = []
        for i in range(m+n):
            thisinterval.append(bound[0]+(i*(bound[1]-bound[0])/(m+n)))
        thisinterval.append(bound[1])
        intervals.append(thisinterval)

    # creating cube of cells with true and false values
    cells = np.full((m+n, nparams), False)
    for i, row in enumerate(cells):     # loop over rows
        for j, cell in enumerate(row):  # loop over columns
            for point in temppoints:    # loop over points to match (shifted by 1 to skip index)
                if point[j+1] >= intervals[j][i] and point[j+1] < intervals[j][i+1]:
                    cells[i][j] = True

    # getting list of cells that are still empty
    emptyCellIDs = []
    for i in range(len(cells[0])):      # loop over parameters
        tempCellIDs = []
        for j in range(len(cells)):     # loop over cells for each parameter
            if cells[j][i] == False:
                tempCellIDs.append(j)
        emptyCellIDs.append(tempCellIDs)        

    # creating points to fill holes
    newpoints = []
    index = 0
    for i in range(n):              # loop for every point to add
        newparams = [index]
        for j in range(nparams):    # loop over each parameter
            # getting a random cell ID to pull parameter value from
            r = random.randint(0,len(emptyCellIDs[j])-1)
            cellID = emptyCellIDs[j][r]
            emptyCellIDs[j].remove(cellID)

            # picking randomly in that cells range
            newparam = (random.random()*(intervals[j][cellID+1]-intervals[j][cellID]))+intervals[j][cellID]
            newparams.append(newparam)

        newpoints.append(newparams)
        index = index+1

    return newpoints

# write the the design out with the header provided
def writeDesign(points, ECM, filename):
    header = "#"
    keysList = [i for s in [d.keys() for d in points[0]] for i in s]
    keysList = keysList[1:] # removing index column
    for key in keysList:
        header = header + " " + key

    newqvirlines = ["# "+ECM+"\n","# Version 1.0\n",header+"\n"]
    for point in points:
        line = ""
        for i, param in enumerate(point[1:]):
            line = line + str(param.get(keysList[i])) + " "
    
        newqvirlines.append(line+"\n")

    qvirCode = open(filename,'w')
    qvirCode.writelines(newqvirlines)
    qvirCode.close()

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