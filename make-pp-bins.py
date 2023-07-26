#makes set of bayesian design points for pp collisions
#first argument is Ecm, with 2760 being treated as a primary run

import os
import multiprocessing as mp
from pickle import TRUE
from pyDOE import lhs 
from datetime import date
from functions import *
from operator import itemgetter
import time

##option reading
system = "LHC"
RHIC = False
Drun = False
Lrun = False
Erun = False
reading = False
rerunning = False
appending = False
design = []
startdir = ""
for i, option in enumerate(sys.argv):
    if "RHIC" in option:
        system = "RHIC"
        RHIC = True
    if "Drun" in option:
        system = "Drun"
        Drun = True
    if "Lrun" in option:
        system = "Lrun"
        Lrun = True
    if "Erun" in option:
        system = "Erun"
        Erun = True
    if "-d" in option:
        reading = True
        design = readDesign(sys.argv[i+1])
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
totaldir = "/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/runs/" + system + "-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"
if rerunning: 
    totaldir = startdir
    design = readDesign(totaldir+"QVir_Analysis/parameters.txt")
    print("Using: "+totaldir+"QVir_Analysis/parameters.txt")


def makexml(bound, parameters, baseDir, xmltemplate, ECM):
    events  = 10000
    #extra events in first bin
    if RHIC:
        events = events*6
    if Drun:
        events = 400000
    if Lrun:
        events = 100000
    if Erun:
        events = 100000
    if bound[0] == 0:
        events = events*2
    if (Lrun or Drun or Erun) and bound[0] <= 10:
        events = events*2

    #appending option and filename
    ogfile = baseDir + "/dat/PP_Bin" + str(bound[0])+ "_" + str(bound[1])
    if appending:
        filename = ogfile+"temp"
    else:
        filename = ogfile

    #building lines for xml file
    eventLine = "  <nEvents> " + str(events) + " </nEvents>"
    lowerLine = "      <pTHatMin>" + str(bound[0]) + "</pTHatMin>\n"
    upperLine = "      <pTHatMax>" + str(bound[1]) + "</pTHatMax>\n"
    fileLine = "  <outputFilename>" + filename + "</outputFilename>\n"
    qLine = "      <Q0>" + str(parameters[1]) + "</Q0>\n"
    virLine = "      <vir_factor>" + str(parameters[2]) + "</vir_factor>\n"
    lambdaline = "    <lambdaQCD> " + str(parameters[3]) + " </lambdaQCD>"
    pionLine = "    <pionWidthScale>" + str(parameters[4]) + "</pionWidthScale>\n"
    kaonLine = "    <kaonWidthScale>" + str(parameters[5]) + "</kaonWidthScale>\n"
    protonLine = "    <protonWidthScale>" + str(parameters[6]) + "</protonWidthScale>\n"
    StoUDline = "      StringFlav:probStoUD = " + str(parameters[7]) + "\n"
    QQtoQline = "      StringFlav:probQQtoQ = " + str(parameters[8]) + "\n"
    pt0refline = "         MultipartonInteractions:pT0Ref = " + str(parameters[9]) +"\n"


    newlines = []
    for line in xmltemplate:
        if "<nEvents>" in line: newlines.append(eventLine)
        elif "<pTHatMin>" in line: newlines.append(lowerLine)
        elif "<pTHatMax>" in line: newlines.append(upperLine)
        elif "<outputFilename>" in line: newlines.append(fileLine)
        elif "<Q0>" in line: newlines.append(qLine)
        elif "<vir_factor>" in line: newlines.append(virLine)
        elif "<lambdaQCD>" in line: newlines.append(lambdaline)
        elif "<pionWidthScale>" in line: newlines.append(pionLine)
        elif "<kaonWidthScale>" in line: newlines.append(kaonLine)
        elif "<protonWidthScale>" in line: newlines.append(protonLine)
        elif "StoUD" in line: newlines.append(StoUDline)
        elif "QQtoQ" in line: newlines.append(QQtoQline)
        elif "pT0Ref" in line: newlines.append(pt0refline)
        else: newlines.append(line)

    xmlname = baseDir+"/xml/pp_"+ECM+"_"+str(parameters[0])+"_"+str(bound[0])+".xml"
    xml = open(xmlname,'w')
    xml.writelines(newlines)
    xml.close()

    output = {'xml': xmlname, 'dat': ogfile}
    return output

#making directory for the runs
def makeDir(index):
    #setting names
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

####################################################################################################################
##Start of main method

##Creating total directory for use
try:
    os.makedirs(totaldir)
except:
    pass

try:
    os.makedirs(totaldir + "QVir_Analysis/")
except:
    pass

##pTHat bounds and xml based on RHIC vs LHC
intervals = []
ECM = ""
if RHIC:
    intervals = [0, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70]
    ECM = "200"
elif Drun or Erun:
    intervals = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100]#, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 230, 250, 270, 290, 310, 330, 350, 400, 450, 500, 550, 600, 1000]
    ECM = "5020"
elif Lrun:
    intervals = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100]#, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 230, 250, 270, 290, 310, 330, 350, 400, 450, 500, 550, 600, 1000]
    ECM = "13000"
else:
    intervals = [0, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 230, 250, 270, 290, 310, 330, 350, 400, 450, 500, 550, 600, 1000]
    ECM = "2760"

xmlname = "/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/config/jetscape_user_pp"+system+".xml"

xmltemplate = open(xmlname,"r")
xmllines = xmltemplate.readlines()

##design points
nsamples = 48
if not reading: design = createDesign(nsamples)

#pthat bounds formed into pairs
pTHatBounds = []
for i in range(len(intervals)-1):
    if i is 0 and intervals[i] == 0:
        pTHatBounds.append((intervals[i],-1))
    else:
        pTHatBounds.append((intervals[i],intervals[i+1]))

##Changing to build directory
os.chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build")

pool = mp.Pool(len(pTHatBounds))

##Loop over design points and make parameter file
newqvirlines = ["# "+ECM+"\n","# Version 1.0\n", "# Parameter q0 virFac lambdaQCD pionScale kaonScale protonScale StoUD QQtoQ pT0Ref\n"]
for index, Q, vir, lqcd, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ, pt0ref in design:
    newqvirlines.append(str(Q)+" "+str(vir)+" "+str(lqcd)+" "+str(pionWidth)+" "+str(kaonWidth)+" "+str(protonWidth)+" "+str(StoUD)+" "+str(QQtoQ)+" "+str(pt0ref)+"\n")

qvirCode = open(totaldir + 'QVir_Analysis/parameters.txt','w')
qvirCode.writelines(newqvirlines)
qvirCode.close()

for index, Q, vir, lQCD, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ, pt0ref in design:
    parameters = [index, Q, vir, lQCD, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ, pt0ref]
    print(parameters)

    ##Making xmls
    baseDir = makeDir(parameters[0])
    output = pool.starmap(makexml, [(bound, parameters, baseDir, xmllines, ECM) for bound in pTHatBounds]) #star map to each set of bounds

    ##Running jetscape for them
    xmls, dats = zip(*map(itemgetter('xml', 'dat'), output)) 
    pool.map(runxml, xmls)

    ##Handling appending runs
    if appending:
        for dat in dats:
            appendcommand = "cat " + dat + "temp.dat >> " + dat + ".dat"
            os.system(appendcommand)
            os.system('rm '+dat+'temp.dat')

pool.close()