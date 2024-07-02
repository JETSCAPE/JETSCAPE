from functions import *
import sys

# making design points
nsamples = sys.argv[1]
design = createPandaDesign(nsamples)
filename = 'designs/' + sys.argv[2] + '.txt'
design.to_csv(filename,index=False)