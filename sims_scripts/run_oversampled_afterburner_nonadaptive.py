import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

print("### Starting afterburner oversampling routine ###")
start_time = time.time()

num_samples = int(sys.argv[1])
print("number of samples : " + str(num_samples))

#number of cores reading the same freezeout surface
num_cores = int(sys.argv[2])
print("Cores available : " + str(num_cores) )

def spawn_afterburner(sample):
    #print('{}: hello from {}'.format( dt.datetime.now(), current_process().name) )
    sample_dir = "sample_" + str(sample)
    os.system( 'mkdir ' + sample_dir )
    os.chdir( sample_dir )
    #link necessary input files to current working dir
    os.system( 'mkdir input' )
    os.chdir( 'input' )
    os.system( 'ln -s ../../surface.dat surface.dat' )
    os.chdir( ".." )

    os.system( 'ln -s ../jetscape_init.xml jetscape_init.xml' )    
    os.system( 'ln -s ../smash_input smash_input' )
    os.system( 'ln -s ../iS3D_parameters.dat iS3D_parameters.dat' )
    os.system( 'ln -s ../PDG PDG' )
    os.system( 'ln -s ../tables tables' )
    os.system( 'ln -s ../deltaf_coefficients deltaf_coefficients' )
    
    non_empty_surf = os.stat("input/surface.dat").st_size
    if (non_empty_surf):
        #run the sampler and afterburner executable, save each stdout to unique file
        os.system( 'SAMPLER_AFTERBURNER > stdout_SAMPLER_AFTERBURNER.txt' )
    else:
        print("***Freezeout surface is empty for this event...***")

    #return to parent dir 
    os.chdir( ".." )

num_launches = num_samples / num_cores

for launch in range(0, num_launches):
    if __name__ == '__main__':
        worker_count = num_cores
        worker_pool = []
        for core in range(worker_count):
            sample = launch * num_cores + core
            p = Process( target = spawn_afterburner, args = (sample,) )
            p.start()
            worker_pool.append(p)
        for p in worker_pool:
            p.join()

print("Oversampling routine finished in " + str( time.time() - start_time) + " sec")
print("Goodbye!")
