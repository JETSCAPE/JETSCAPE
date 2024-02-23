from functions import *
import os
import sys

datadir = "/scratch/user/cameron.parker/projects/STAT/input/data/pp/"
finaldir = "/scratch/user/cameron.parker/projects/STAT/input/data/pp-copy/"
inputdir = sys.argv[1] + "QVir_Analysis/"

#bounds to remove
lower = 6
upper = 10

#indices for range to slice out
i1 = 0
i2 = 0

datafiles = ["2760chargedhads.dat",
             "2760pions.dat",
             "2760kaons.dat",
             "2760protons.dat"]

predictfiles = ["HadronSpectraPrediction.dat",
                "PionSpectraPrediction.dat",
                "KaonSpectraPrediction.dat",
                "ProtonSpectraPrediction.dat"]

for i in range(len(predictfiles)):
    datfile = open(datadir+datafiles[i],"r")
    datlines = datfile.readlines()

    infile = open(inputdir+predictfiles[i],"r")
    inlines = infile.readlines()

    for iline,line in enumerate(datlines):
        if("#" in line): continue

        split = line.split()
        if(float(split[0]) == lower):
            i1 = len(datlines)-iline
            
        if(float(split[1]) <= upper):
            i2 = len(datlines)-iline

    newdatlines = []
    for iline,line in enumerate(datlines):
        if(len(datlines) - iline <= i1 and len(datlines) - iline >= i2):
            continue
        else:
            newdatlines.append(line)

    newdatfile = open(finaldir+datafiles[i],"w")
    newdatfile.writelines(newdatlines)

    newpredlines = []
    for iline,line in enumerate(inlines):
        if(len(inlines) - iline <= i1 and len(inlines) - iline >= i2):
            continue
        else:
            newpredlines.append(line)

    newpredfile = open(inputdir+predictfiles[i]+"2","w")
    newpredfile.writelines(newpredlines)
