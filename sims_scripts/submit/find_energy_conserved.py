import sys
import os
import re
import numpy as np

#import matplotlib.pyplot as plt

dir = sys.argv[1]

part = "Total Energy changed by "

#radii = []
energy_file = open('energy.dat', 'w')

for subdir, dirs, files in os.walk(dir):
    for file in files:
        #print os.path.join(subdir, file)
        filepath = subdir + os.sep + file
        #print(filepath)
        with open(filepath, 'r') as file:
            file_txt = file.read()
            #print(file_txt)
            before_keyword, keyword, after_keyword = file_txt.partition(part)
            if (len(after_keyword) > 0):
                radius = after_keyword.split()[0]
                #radii.append(float(radius))
                rad_file.write(str(radius) + '\n')

#radii = np.array(radii)
#np.savetxt('radii.dat', radii)
rad_file.close()

#plt.hist(radii)
#plt.savefig('radii.png')
