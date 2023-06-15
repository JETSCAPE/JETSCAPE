import os
import sys

#gets list of directories in the target directory
def getDats(dir):
    startdir = os.getcwd()
    print(dir)
    os.chdir(dir+'/dat')
    files = os.listdir(".")
    
    os.chdir(startdir)
    
    files.sort()  #to make sure the directories are sorted in case it doesnt complete
    return files

def getDirs(dir):
    startdir = os.getcwd()
    os.chdir(dir)
    files = os.listdir(".")

    directories = []
    for item in files:
        if os.path.isdir(item) and item != 'QVir_Analysis': directories.append(item)

    os.chdir(startdir)

    
    directories.sort()  #to make sure the directories are sorted in case it doesnt complete
    return directories

#setting directory for analysis
basedir = sys.argv[1]
dirs = getDirs(basedir)

for analysisDir in dirs:
    files = getDats(basedir+analysisDir)

    for file in files:
        datfile = open(basedir+analysisDir+'/dat/'+file,'r')
        datlines = datfile.readlines()
        eventline = 'no events'
        eventcount = 0

        for line in datlines:
            if 'Event' in line:
                eventline = line
                eventcount = eventcount + 1

        print(str(file) + ": " + str(eventcount))
        