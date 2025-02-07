from functions import *
import sys

# making design points
nsamples = int(sys.argv[1])
design = createPandaDesign(nsamples)
filename = 'designs/' + sys.argv[2] + '.txt'
design.to_csv(filename,index=False)