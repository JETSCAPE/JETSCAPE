#running over points for pp collisions

import os
import multiprocessing as mp
from functions import *

totdir = sys.argv[1]
dirs = getDirs(totdir)

# selecting range to run over
start = 0
finish = len(dirs)
dirs = intSort(dirs)

if len(sys.argv) == 4:
    start = int(sys.argv[2])
    finish = int(sys.argv[3]) + 1

print("From " + str(start) + " to " + str(finish) + " in " + totdir + ":")

# Changing to build directory
os.chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build")
pool = mp.Pool(48)

# looping over points
for dir in dirs[start:finish]:
    print("Running " + dir + ":")
    pointdir = totdir+"points/"+dir
    unzipxmls(pointdir)
    xmls = readXmls(pointdir)
    pool.map(runxml, xmls)
    zipxmls(pointdir)

    # concatonating all soft bins together
    if "_0_" in xmls[0] and "_0_" in xmls[1]:
        dats = concatDats(totdir+dir)

pool.close()
    