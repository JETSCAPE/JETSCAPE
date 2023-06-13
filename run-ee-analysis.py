#runs analysis over ee bayesian set of design points. takes target dir as an argument

import os
import sys
from functions import *

#setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)

os.chdir("/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/build/")

for directory in directories:
    baseDir = "/scratch/user/cameron.parker/JETSCAPE-COMP-HH_colorrecomb/" + analysisDir + directory
    cmd = "./ee-analysis-spectra " + baseDir
    update(cmd)

    os.system(cmd)

os.system("./ee-comparison ../" + analysisDir)