import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

print("### Starting multithread sampling routine ###")
start_time = time.time()

#number of samples
nsamples = int(sys.argv[1])
print("Number of samples : " + str(nsamples) )
#number of cores reading the same freezeout surface
ncores = int(sys.argv[2])
print("Cores available : " + str(ncores) )

def spawn_sampler(sample):
    sample_dir = "sample_" + str(sample)
    os.system( 'mkdir ' + sample_dir )
    os.chdir( sample_dir )
    #link necessary input files to current working dir
    os.system( 'ln -s ../jetscape_init.xml jetscape_init.xml' )
    os.system( 'ln -s ../surface.dat surface.dat' )
    os.system( 'ln -s ../iSS_parameters.dat iSS_parameters.dat' )
    os.system( 'ln -s ../music_input music_input' )
    os.system( 'ln -s ../iSS_tables iSS_tables' )
    os.system( 'ln -s ../iSSSamplingOnly iSSSamplingOnly' )
    #run the sampler
    os.system( './iSSSamplingOnly' )
    #return to parent dir 
    os.chdir( ".." )

num_launches = (nsamples / ncores) + 1

#spawn the jobs
for launch in range(0, num_launches):
    if __name__ == '__main__':
        worker_count = ncores
        worker_pool = []
        for core in range(worker_count):
            sample = launch * ncores + core
            p = Process( target = spawn_sampler, args = (sample,) )
            p.start()
            worker_pool.append(p)
        for p in worker_pool:
            p.join()

print("Oversampling routine finished in " + str( time.time() - start_time) + " sec")
print("Goodbye!")
