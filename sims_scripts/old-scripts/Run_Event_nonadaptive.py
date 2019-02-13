import csv
import os
import sys
from multiprocessing import Process, current_process
import datetime as dt
import time

print( "### Starting Event ###" )
start_time = time.time()
num_samples = int(sys.argv[1])

print( "### Number of samples : " + str(num_samples) + " ###" )

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

#first call the executable which runs TRENTO + FREESTREAMING + HYDRO
#this writes surface.dat to cwd
os.system( 'TRENTO_FS_HYDRO > stdout_TRENTO_FS_HYDRO.txt' )

#check if freezeout surface is empty
#if so, skip sampler and afterburner
non_empty_surf = os.stat("surface.dat").st_size
if (non_empty_surf):
    #call script that oversamples surface and runs afterburners
    os.system( 'python run_oversampled_afterburner_nonadapative.py ' + str(num_samples) )
else:
    print("***Freezeout surface is empty for this event...***")

end_time = time.time()

print( "Event finished in " + str( end_time - start_time) + " sec" )
