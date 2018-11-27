import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

print("### Starting afterburner oversampling routine ###")
start_time = time.time()

# stop spawning sampling/afterburner events when
# the total number of hadrons summed over all samples
# is greater than some lower bound
min_num_particles = int(sys.argv[1])
print("Minimum number of total particles : " + str(min_num_particles))

#number particles sampled in first set
num_particles_sampled = 0

#number of cores reading the same freezeout surface
num_cores = int(sys.argv[2])
print("Cores available : " + str(num_cores) )

def spawn_afterburner(sample):
    #print('{}: hello from {}'.format( dt.datetime.now(), current_process().name) )
    sample_dir = "sample_" + str(sample)
    os.system( 'mkdir ' + sample_dir )
    os.chdir( sample_dir )
    #link necessary input files to current working dir
    os.system( 'ln -s ../jetscape_init.xml jetscape_init.xml' )
    os.system( 'ln -s ../surface.dat surface.dat' )
    os.system( 'ln -s ../music_input music_input' )
    os.system( 'ln -s ../smash_input smash_input' )
    os.system( 'ln -s ../iSS_parameters.dat iSS_parameters.dat' )
    os.system( 'ln -s ../iSS_tables iSS_tables' )
    os.system( 'ln -s ../SAMPLER_AFTERBURNER SAMPLER_AFTERBURNER' )
    #run the sampler and afterburner executable
    os.system( './SAMPLER_AFTERBURNER' )
    #return to parent dir 
    os.chdir( ".." )

def get_number_particles(sample):
    sample_dir = "sample_" + str(sample)
    #particle_list_file = open(sample_dir + '/final_smash_hadrons.dat', 'rb')
    #reader = csv.reader(particle_list_file)
    #first header line is number of particles in list
    #num_particles = reader.next()
    #return int(num_particles[0])

    #this method assumes fixed length of header
    with open(sample_dir + '/final_smash_hadrons.dat', 'rb') as f:
        line_count = 0
        for line in f:
            line_count += 1
    return line_count - 1


#spawn the first set of jobs
if __name__ == '__main__':
    worker_count = num_cores
    worker_pool = []
    for sample in range(worker_count):
        p = Process( target = spawn_afterburner, args = (sample,) )
        p.start()
        worker_pool.append(p)
    for p in worker_pool:
        p.join()

#get the number of particles produced in the first set of samples
for sample in range(0, num_cores):
    num_particles = get_number_particles(sample)
    num_particles_sampled += num_particles

print("Number of particles sampled in first launch : " + str(num_particles_sampled) )
avg_num_per_sample = num_particles_sampled / num_cores
print("Average number of particles per sample : " + str(avg_num_per_sample) )

#get the total number of events we need based on the avg
num_jobs = min_num_particles / avg_num_per_sample + 1

print("Total number of jobs necessary : " + str(num_jobs) )
num_launches = ( num_jobs - 1 ) / num_cores
print("Number of launches necessary to meet minimum : " + str(num_launches) )

for launch in range(1, num_launches + 1):
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
