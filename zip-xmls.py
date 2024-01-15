from functions import *

totdir = sys.argv[1]
dirs = getDirs(totdir)
dirs = intSort(dirs)

# looping over points
for dir in dirs:
    zipxmls(totdir+"points/"+dir)