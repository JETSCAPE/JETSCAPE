import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

num_events = int(sys.argv[1])
print("Number of hydro events : " + str(num_events))

num_cores = int(sys.argv[2])
print("Number of parallel cores : " + str(num_cores))

num_samples = int(sys.argv[3])
print("number of samples per surface : " + str(num_samples))

print("### Starting multiple hydro jobs on node ###")
start_time = time.time()


def spawn_event(event):
    event_dir = "event_" + str(event)
    os.system( 'mkdir ' + event_dir )
    os.chdir( event_dir )
    #link necessary input files to current working dir
    #files for TRENTO+FS+HYDRO
    os.system( 'ln -s ../jetscape_init.xml jetscape_init.xml' )
    os.system( 'ln -s ../freestream_input freestream_input' )
    os.system( 'ln -s ../music_input music_input' )
    os.system( 'ln -s ../EOS EOS' )
    #files for SAMPLER+AFTERBURNER
    os.system( 'ln -s ../run_oversampled_afterburner_nonadaptive.py run_oversampled_afterburner_nonadaptive.py' )
    os.system( 'ln -s ../smash_input smash_input' )
    os.system( 'ln -s ../iS3D_parameters.dat iS3D_parameters.dat' )
    os.system( 'ln -s ../PDG PDG' )
    os.system( 'ln -s ../tables tables' )
    os.system( 'ln -s ../deltaf_coefficients deltaf_coefficients' )
    
    #run the TRENTO+FS+HYDRO executable and save stdout to unique file
    os.system( 'TRENTO_FS_HYDRO > stdout_TRENTO_FS_HYDRO.txt' )

    #now run the oversampling and afterburners, only call one core, other cores are busy with other events...
    os.system( 'python run_oversampled_afterburner_nonadaptive.py ' + str(num_samples) + ' 1' )

    #return to parent dir 
    os.chdir( ".." )


#determine number of launches necessary to meet number of events
num_launches = (num_events - 1) / num_cores + 1

for launch in range(0, num_launches):
    if __name__ == '__main__':
        worker_count = num_cores
        if (launch == num_launches - 1):
            worker_count = num_events - launch * num_cores
        worker_pool = []
        for core in range(worker_count):
            event = launch * num_cores + core
            p = Process( target = spawn_event, args = (event,) )
            p.start()
            worker_pool.append(p)
        for p in worker_pool:
            p.join()


print("Events multicore routine finished in " + str( time.time() - start_time) + " sec")
print("Goodbye!")
