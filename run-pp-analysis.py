#runs the analysis over an entire set of design points. takes target dir and either "full" or "smoothing" as inputs

import os
import sys
import multiprocessing as mp
from functions import *

#setting run options
ECM = "2760"
RHICrun = False
LHC900 = False
LHC7000 = False
smoothrun = False
setStart = False
parallel = False
startdir = 0
for i, option in enumerate(sys.argv):
    if "RHIC" in option:
        RHICrun = True
        ECM = "200"
    if "LHC900" in option:
        LHC900 = True
        ECM = "900"
    if "LHC7000" in option:
        LHC7000 = True
        ECM = "7000"
    if "smooth" in option:
        smoothrun = True
    if "-s" in option and option.startswith("-"):
        setStart = True
        startdir = sys.argv[i+1]
    if "-p" in option and option.startswith("-"):
        parallel = True

#setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
intSort(directories)

#skipping ones we already did if start specified
if setStart:
    newdirectories = []
    skip = True
    for directory in directories:
        if startdir == directory:
            skip = False
        if skip == False:
            newdirectories.append(directory)
    directories = newdirectories


os.chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build/")

def run(directory):
    if analysisDir.startswith("/"):
        baseDir = analysisDir +  "points/" + directory
    else:
        baseDir = "/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/" + analysisDir + "points/" + directory
    
    if smoothrun:
        cmd = "./pp-"+ECM+"-smoothing " + baseDir
    else:
        cmd = "./pp-"+ECM+"-spectra " + baseDir

    update(cmd)
    os.system(cmd)
    print("Point finished.")
    return

#Directory loop
if parallel:
    pool = mp.Pool(len(directories))
    pool.map(run,directories)
    pool.close()
else:
    for directory in directories:
        run(directory)

#Total analysis
os.system("./pp-"+ECM+"-comparison ../" + analysisDir)