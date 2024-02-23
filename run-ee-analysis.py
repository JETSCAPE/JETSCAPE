#runs analysis over ee bayesian set of design points. takes target dir as an argument

import os
import sys
from functions import *
import multiprocessing as mp

#option reading
parallel = False
for i, option in enumerate(sys.argv):
    if "-p" in option and option.startswith("-"):
        parallel = True

#setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)

os.chdir("/scratch/user/cameron.parker/newJETSCAPE/build/")

def run(directory):
    if "/scratch" in analysisDir:
        baseDir = analysisDir + "points/" + directory
    else:
        baseDir = "/scratch/user/cameron.parker/newJETSCAPE/" + analysisDir + "points/" + directory
    cmd = "./ee-analysis-spectra " + baseDir
    update(cmd)

    os.system(cmd)

#Directory loop
if parallel:
    pool = mp.Pool(48)
    pool.map(run,directories)
    pool.close()
else:
    for directory in directories:
        run(directory)

os.system("./ee-comparison ../" + analysisDir)