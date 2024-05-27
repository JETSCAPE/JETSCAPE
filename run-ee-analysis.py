#runs analysis over ee bayesian set of design points. takes target dir as an argument

import os
import sys
from functions import *
import multiprocessing as mp
import htcondor
import classad

# option reading
parallel = False
for i, option in enumerate(sys.argv):
    if "-p" in option and option.startswith("-"):
        parallel = True

# job initialization
condorjob = analysisjob
schedd = htcondor.Schedd()                   # get the Python representation of the scheduler

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)

os.chdir("/data/rjfgroup/rjf01/cameron.parker/builds/JETSCAPE/build/")

def run(directory):
    if "/scratch" in analysisDir:
        baseDir = analysisDir + "points/" + directory
    else:
        baseDir = "/data/rjfgroup/rjf01/cameron.parker/runs/" + analysisDir + "points/" + directory
    cmd = "./ee-analysis-spectra " + baseDir
    update(cmd)

    os.system(cmd)

# Directory loop
if parallel:
    dirinput = [{"dir": analysisDir+"points/"+dir} for dir in directories]
    condorjob["batch_name"] = analysisDir.split('/')[-2]+"-analysis"
    submit_result = schedd.submit(condorjob, itemdata = iter(dirinput))
    print(dirinput)
else:
    for directory in directories:
        run(directory)