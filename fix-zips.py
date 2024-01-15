from functions import *
import shutil
import os

totdir = sys.argv[1]
dirs = getDirs(totdir)
dirs = intSort(dirs)

# looping over points
for dir in dirs[59:]:
    pointdir = totdir+"points/"+dir
    print("Fixing "+pointdir)
    sourcedir = pointdir+"/xml/runs/LHC-1-2/points/"+dir+"/xml"
    destination = pointdir+"/xml"
    files = os.listdir(sourcedir)

    for file in files:
        src_path = os.path.join(sourcedir, file)
        dst_path = os.path.join(destination, file)
        shutil.move(src_path, dst_path)

    shutil.rmtree(pointdir+"/xml/runs")
