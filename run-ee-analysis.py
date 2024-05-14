#runs analysis over ee bayesian set of design points. takes target dir as an argument

import os
import sys
from functions import *
import multiprocessing as mp
import htcondor
import classad

#option reading
parallel = False
for i, option in enumerate(sys.argv):
    if "-p" in option and option.startswith("-"):
        parallel = True

#job initialization
testjob = htcondor.Submit({
    "executable": "/data/rjfgroup/rjf01/cameron.parker/builds/JETSCAPE/build/ee-analysis-spectra",
    "arguments": "$(dir)",          # we will pass in the value for this macro via itemdata
    "output": "/data/rjfgroup/rjf01/cameron.parker/condor/cat-$(ProcId).out",
    "error": "/data/rjfgroup/rjf01/cameron.parker/condor/cat-$(ProcId).err",
    "log": "/data/rjfgroup/rjf01/cameron.parker/condor/cat.log",
    "request_cpus": "1",
    "request_memory": "500MB",
    "request_disk": "500MB",
})
schedd = htcondor.Schedd()                   # get the Python representation of the scheduler

#setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)

os.chdir("/scratch/user/cameron.parker/projects/JETSCAPE/build/")

def run(directory):
    if "/scratch" in analysisDir:
        baseDir = analysisDir + "points/" + directory
    else:
        baseDir = "/scratch/user/cameron.parker/projects/JETSCAPE/" + analysisDir + "points/" + directory
    cmd = "./ee-analysis-spectra " + baseDir
    update(cmd)

    os.system(cmd)

#Directory loop
if parallel:
    dirinput = [{"dir": dir} for dir in directories]
    submit_result = schedd.submit(testjob, itemdata = iter(dirinput))
else:
    for directory in directories:
        run(directory)