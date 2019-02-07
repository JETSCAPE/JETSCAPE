import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

print( "### Starting Event ###" )
start_time = time.time()

# stop spawning sampling/afterburner events when the total number 
# of hadrons summed over all samples is greater than some lower bound
min_num_particles = int(sys.argv[1])

print( "### Minimum number Particles : " + str(min_num_particles) + " ###" )

#number of cores reading the same freezeout surface
num_cores = int(sys.argv[2])

print( "### Cores available : " + str(num_cores) + " ###" )

#first call the executable which runs TRENTO + FREESTREAMING + HYDRO
#this writes surface.dat to cwd
os.system( 'TRENTO_FS_HYDRO' )

#now call the afterburner oversampling script
os.system( 'python run_oversampled_afterburner.py ' + str(min_num_particles) + ' ' + str(num_cores) )

end_time = time.time()

print( "Event finished in " + str( end_time - start_time) + " sec" )
