import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

print("### Starting multiple hydro jobs on node ###")
start_time = time.time()

num_cores = int(sys.argv[1])
print("Number of parallel instances of hydro : " + str(num_cores))

def spawn_hydro(event):
    event_dir = "event_" + str(event)
    os.system( 'mkdir ' + event_dir )
    os.chdir( event_dir )
    #link necessary input files to current working dir
    os.system( 'ln -s ../jetscape_init.xml jetscape_init.xml' )
    os.system( 'ln -s ../freestream_input freestream_input' )
    os.system( 'ln -s ../music_input music_input' )
    os.system( 'ln -s ../EOS EOS' )
    #having these will be useful later for sampling and afterburner
    os.system( 'ln -s ../smash_input smash_input' )
    os.system( 'ln -s ../iS3D_parameters.dat iS3D_parameters.dat' )
    os.system( 'ln -s ../PDG PDG' )
    os.system( 'ln -s ../tables tables' )
    os.system( 'ln -s ../deltaf_coefficients deltaf_coefficients' )
    
    #run the TRENTO+FS+HYDRO executable
    os.system( 'TRENTO_FS_HYDRO' )
    #return to parent dir 
    os.chdir( ".." )

#spawn the first set of jobs
if __name__ == '__main__':
    worker_count = num_cores
    worker_pool = []
    for event in range(worker_count):
        p = Process( target = spawn_hydro, args = (event,) )
        p.start()
        worker_pool.append(p)
    for p in worker_pool:
        p.join()

print("Hydro multicore routine finished in " + str( time.time() - start_time) + " sec")
print("Goodbye!")
