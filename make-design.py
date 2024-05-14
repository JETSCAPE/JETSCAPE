from functions import *

# making design points
nsamples = 500
design = createPandaDesign(nsamples)
design.to_csv('designs/totaldesign.txt',index=False)