import os
import sys
from pyDOE import lhs 
import numpy as np
import random

#gets list of directories in the target directory
def getDirs(dir):
    startdir = os.getcwd()
    os.chdir(dir)
    files = os.listdir(".")

    directories = []
    for item in files:
        if os.path.isdir(item) and item != 'QVir_Analysis': directories.append(item)

    os.chdir(startdir)

    
    directories.sort()  #to make sure the directories are sorted in case it doesnt complete
    return directories

#to see what ran last if runs dont complete
def update(cmd):
    outfile = open('last-ran.txt','a')
    outfile.write(cmd+'\n')
    outfile.close()

    print(cmd)

#reading design file from an existing run
def readDesign(designfile):
    design = []

    #opening file
    file = open(designfile, 'r')
    lines = file.readlines()

    index = 0
    for line in lines:
        #skipping comments
        if '#' in line:
            continue
        
        #splitting line and converting to floats
        params = line.split()
        newparams = [index]
        for param in params:
            newparams.append(float(param))
        design.append(newparams)
        index = index + 1

    return design

##ranges of form [(start,end),(start,end)]
def createDesign(points):
    ##initialization
    rawdesign = lhs(8, samples = points)
    index = 0
    design = []
    print(rawdesign)

    #scaling to param ranges
    for Q, vir, lQCD, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ in rawdesign:
        ##Rescaling LHS parameters to appropriate ranges
        Qnew = (Q*2.1 + 0.9)
        virnew = (vir*.9 + .1)
        lQCDnew = (lQCD*0.3 + 0.1)
        pionWidthnew = (pionWidth*1.5 + 0.5)
        kaonWidthnew = (kaonWidth*1.5 + 0.5)
        protonWidthnew = (protonWidth*1.5 + 0.5)
        StoUDnew = (StoUD*0.3+0.2)
        QQtoQnew = (QQtoQ*0.13+0.07)
        design.append([index, Qnew, virnew, lQCDnew, pionWidthnew, kaonWidthnew, protonWidthnew, StoUDnew, QQtoQnew])
        index = index + 1

    return design

##augment latin hypercube
def augmentDesign(points, bounds, n):
    m = len(points)
    nparams = len(points[0])-1
    temppoints = points

    #setting up intervals for hypercube cells
    intervals = []
    for bound in bounds:    #loop over parameters
        thisinterval = []
        for i in range(m+n):
            thisinterval.append(bound[0]+(i*(bound[1]-bound[0])/(m+n)))
        thisinterval.append(bound[1])
        intervals.append(thisinterval)

    #creating cube of cells with true and false values
    cells = np.full((m+n, nparams), False)
    for i, row in enumerate(cells):     #loop over rows
        for j, cell in enumerate(row):  #loop over columns
            for point in temppoints:    #loop over points to match (shifted by 1 to skip index)
                if point[j+1] >= intervals[j][i] and point[j+1] < intervals[j][i+1]:
                    cells[i][j] = True

    #getting list of cells that are still empty
    emptyCellIDs = []
    for i in range(len(cells[0])):      #loop over parameters
        tempCellIDs = []
        for j in range(len(cells)):     #loop over cells for each parameter
            if cells[j][i] == False:
                tempCellIDs.append(j)
        emptyCellIDs.append(tempCellIDs)        

    #creating points to fill holes
    newpoints = []
    index = 0
    for i in range(n):              #loop for every point to add
        newparams = [index]
        for j in range(nparams):    #loop over each parameter
            #getting a random cell ID to pull parameter value from
            r = random.randint(0,len(emptyCellIDs[j])-1)
            cellID = emptyCellIDs[j][r]
            emptyCellIDs[j].remove(cellID)

            #picking randomly in that cells range
            newparam = (random.random()*(intervals[j][cellID+1]-intervals[j][cellID]))+intervals[j][cellID]
            newparams.append(newparam)

        newpoints.append(newparams)
        index = index+1

    return newpoints

#write the the design out with the header provided
def writeDesign(points, header, filename):
    newqvirlines = [header]
    for index, Q, vir, lqcd, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ in points:
        newqvirlines.append(str(Q)+" "+str(vir)+" "+str(lqcd)+" "+str(pionWidth)+" "+str(kaonWidth)+" "+str(protonWidth)+" "+str(StoUD)+" "+str(QQtoQ)+"\n")

    qvirCode = open(filename,'w')
    qvirCode.writelines(newqvirlines)
    qvirCode.close()

#runs the xml in jetscape
def runxml(xmlname):
    runjetscapecmd = "./runJetscape " + xmlname
    update(runjetscapecmd)
    os.system(runjetscapecmd)