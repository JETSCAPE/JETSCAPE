#makes set of bayesian design points for ee collisions

import os
import multiprocessing as mp
from pyDOE import * 
from datetime import date
from functions import *

#date and file initiation
todays_date = date.today()
totaldir = "/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/runs/LEP-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"

try:
    os.makedirs(totaldir)
except:
    pass

try:
    os.makedirs(totaldir + "QVir_Analysis/")
except:
    pass

#methods for runnning
def binrun(designpoint, xmltemplate):
    ##making directory for this designpoint
    baseDir = makeDir(designpoint[0])

    ##building lines for xml file
    fileLine = "    <outputFilename>" + baseDir + "/run" + "</outputFilename>\n"
    qLine = "            <Q0>" + str(designpoint[1]) + "</Q0>\n"
    virLine = "            <vir_factor>" + str(designpoint[2]) + "</vir_factor>\n"
    lambdaline = "        <lambdaQCD> " + str(designpoint[3]) + " </lambdaQCD>"
    pionLine = "        <pionWidthScale>" + str(designpoint[4]) + "</pionWidthScale>\n"
    kaonLine = "        <kaonWidthScale>" + str(designpoint[5]) + "</kaonWidthScale>\n"
    protonLine = "        <protonWidthScale>" + str(designpoint[6]) + "</protonWidthScale>\n"
    StoUDline = "        <StoUD>" + str(designpoint[7]) + "</StoUD>\n"
    QQtoQline = "        <QQtoQ>" + str(designpoint[8]) + "</QQtoQ>\n"

    newlines = []
    for line in xmltemplate:
        if "<outputFilename>" in line: newlines.append(fileLine)
        elif "<Q0>" in line: newlines.append(qLine)
        elif "<vir_factor>" in line: newlines.append(virLine)
        elif "<lambdaQCD>" in line: newlines.append(lambdaline)
        elif "<pionWidthScale>" in line: newlines.append(pionLine)
        elif "<kaonWidthScale>" in line: newlines.append(kaonLine)
        elif "<protonWidthScale>" in line: newlines.append(protonLine)
        elif "<StoUD>" in line: newlines.append(StoUDline)
        elif "<QQtoQ>" in line: newlines.append(QQtoQline)
        else: newlines.append(line)

    xmlname = baseDir + "/config.xml"
    xml = open(xmlname,'w')
    xml.writelines(newlines)
    xml.close()

    ##running for the parameters set in the xml file
    runxml(xmlname)

    analysiscmd = "./ee-analysis-spectra " + baseDir
    #os.system(analysiscmd)

##making directory for the runs
def makeDir(index):
    ##setting names
    baseDir = totaldir + str(index)

    ##making dirs
    try:
        os.makedirs(baseDir)
    except:
        pass

    ##returning base directory name
    return baseDir

######################################################################################################################################
##start of main method

##option reading
reading = False
if len(sys.argv) == 2:
    reading = True

##parameters and bounds
points = 48
tempdesign = lhs(8, samples = points)
#tempdesign.sort()
design = []
index = 0
for Q, vir, lQCD, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ in tempdesign:
    ##Rescaling LHS parameters to appropriate ranges
    Qnew = (Q*2.1 + 0.9)
    virnew = -(vir*0.9 + 0.1)
    lQCDnew = (lQCD*0.3 + 0.1)
    pionWidthnew = (pionWidth*1.5 + 0.5)
    kaonWidthnew = (kaonWidth*1.5 + 0.5)
    protonWidthnew = (protonWidth*1.5 + 0.5)
    StoUDnew = (StoUD*0.3+0.2)
    QQtoQnew = (QQtoQ*0.13+0.07)
    design.append([index, Qnew, virnew, lQCDnew, pionWidthnew, kaonWidthnew, protonWidthnew, StoUDnew, QQtoQnew])
    index = index+1

if reading:
    design = readDesign(sys.argv[1])

#reading xml template
xmltemplate = open("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/config/jetscape_user.xml","r")
xmllines = xmltemplate.readlines()

##Changing to build directory
os.chdir("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/build")

##setting parameters to run analysis over
newqvirlines = [ '# Version 1.0\n', '# Parameter q0 virFac lambdaQCD pionScale kaonScale protonScale StoUD QQtoQ\n']
for index, Q, vir, lqcd, pionWidth, kaonWidth, protonWidth, StoUD, QQtoQ in design:
    newqvirlines.append(str(Q)+" "+str(vir)+" "+str(lqcd)+" "+str(pionWidth)+" "+str(kaonWidth)+" "+str(protonWidth)+" "+str(StoUD)+" "+str(QQtoQ)+"\n")

qvirCode = open(totaldir + '/QVir_Analysis/parameters.txt','w')
qvirCode.writelines(newqvirlines)
qvirCode.close()

##Running jetscape for each set of parameters with multiprocessing
pool = mp.Pool(points)
pool.starmap(binrun, [(point, xmllines) for point in design])
pool.close()

##running total analysis
#os.system("./ee-comparison " + totaldir)