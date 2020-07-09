import argparse
from contextlib import contextmanager
import datetime
from itertools import chain, groupby, repeat
import logging
import math
import os
import pickle
import signal
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
import shutil

from eccentricities import calc_ecc
#from find_tau_fs import calc_tau_fs

# fully specify numeric data types, including endianness and size, to
# ensure consistency across all machines
float_t = '<f8'
int_t = '<i8'
complex_t = '<c16'

# species (name, ID) for identified particle observables
species = [
	('pion', 211),
	('kaon', 321),
	('proton', 2212),
	('Lambda', 3122),
	('Sigma0', 3212),
	('Xi', 3312),
	('Omega', 3334),
	('phi', 333),
	#('d', 1000010020),
]
pi_K_p = [
	('pion', 211),
	('kaon', 321),
	('proton', 2212),
]

Qn_species = [
	('pion', 211),
	('kaon', 321),
	('proton', 2212),
	('Sigma', 3222),
	('Xi', 3312),
	#('d', 1000010020),
	('Ch', None)
]

#old fine
# Qn_diff_pT_cuts=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.2,3.4,3.6,3.8,4.]
#updated fine
#Qn_diff_pT_cuts=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.2,3.4,3.6,3.8,4.,10.]
#coarse
Qn_diff_pT_cuts=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]


Qn_diff_NpT = len(Qn_diff_pT_cuts)-1
Nharmonic = 8
Nharmonic_diff = 5

# Particle datatype
particle_dtype = [
                                ('ID', int),
                                ('charge', int),
                                ('pT', float),
                                ('ET', float),
                                ('mT', float),
                                ('phi', float),
                                ('y', float),
                                ('eta', float),
                        ]

# results "array" (one element)
# to be overwritten for each event
#five entries for five 'delta-f' options: 0,1,2,3 run with smash, 4 runs with urqmd
expt_type = 'STAR'
number_of_viscous_corrections=4
result_dtype=[
('initial_entropy', float_t, 1),
('impact_parameter', float_t, 1),
('npart', float_t, 1),
(expt_type,
        [
                ('nsamples', int_t, 1),
                # 1) dNch/deta, eta[-0.5, 0.5], charged
                ('dNch_deta', float_t, 1),
                # 2) dET/deta, eta[-0.6, 0.6]
                ('dET_deta', float_t, 1),
                # 3.1) The Tmunu observables, eta[-0.6, 0.6]
                ('Tmunu', float_t, 10),
                # 3.2) The Tmunu observables, eta[-0.5, 0.5], charged
                ('Tmunu_chg', float_t, 10),
                # 4.1) identified particle yield
                ('dN_dy',       [(name, float_t, 1) for (name,_) in species], 1),
                # 4.2) identified particle <pT>
                ('mean_pT', [(name, float_t, 1) for (name,_) in species], 1),
                # 5.1) pT fluct, pT[0.15, 2], eta[-0.8, 0.8], charged
                ('pT_fluct_chg', [      ('N', int_t, 1),
                                                        ('sum_pT', float_t, 1),
                                                        ('sum_pT2', float_t, 1)], 1),
                # 5.2) pT fluct, pT[0.15, 2], eta[-0.8, 0.8], pi, K, p
                ('pT_fluct_pid', [      (name, [        ('N', int_t, 1),
                                                                                ('sum_pT', float_t, 1),
                                                                                ('sum_pT2', float_t, 1)], 1     )
                                                          for (name,_) in pi_K_p        ], 1),
                # 6) Q vector, pT[0.2, 5.0], eta [-0.8, 0.8], charged
                ('flow', [      ('N', int_t, 1),
                                        ('Qn', complex_t, Nharmonic)], 1),
        ], number_of_viscous_corrections),
# Q vector, diff-flow, identified charged hadrons
('d_flow_pid', [(name, [('N', int_t, Qn_diff_NpT),
                                                ('Qn', complex_t, [Qn_diff_NpT, Nharmonic_diff])], 1)
                                for (name,_) in Qn_species      ], number_of_viscous_corrections),
]

# A function that calculates event-planes given a list of azimuthal angles
def calculate_eventplanes(phi, Nharmonic):
        V_n = {n: float_t for n in range(1,Nharmonic+1)}
        Psi_n =  {n: float_t for n in range(1,Nharmonic+1)}
        N = phi.size
        Npairs = N*(N-1.)
        for n in range(1, Nharmonic+1):
                Qx = np.cos(n*phi).sum()
                Qy = np.sin(n*phi).sum()
                V_n[n] = np.sqrt( (Qx**2 + Qy**2 - N)/pairs )
                Psi_n[n] = np.arctan2(Qy, Qx)/n
        return V_n, Psi_n

# A function that calculates double differential q-vectors
# given a list of pT-bins and a list of eta-bins
# for example:
#        pTbins = np.array([[0,1],[1,2],[2,3]])
#       etabins = np.array([[-3,-1],[-1,1],[1,3]])
def calc_Qvector(pT, phi, eta, pTbins, etabins, Nharmonic):
        shape = [pTbins.shape[0], etabins.shape[0]]
        # calculate Qn
        results = {'N': np.zeros(shape, dtype=int_t),
                           'Qn': np.zeros([*shape, Nharmonic], dtype=complex_t)}
        for i, (pT_l, pT_h) in enumerate(pTbins):
                pT_cut = (pT_l<=pT) & (pT<=pT_h)
                for j, (eta_l, eta_h) in enumerate(etabins):
                        eta_cut = (eta_l<=eta) & (eta<=eta_h)
                        # apply cut
                        sub_phi = phi[pT_cut & eta_cut]
                        results['N'][i,j] = sub_phi.size

                        results['Qn'][i,j] = np.array( [ np.exp(1j*n*sub_phi).sum() for n in range(1, 1+Nharmonic) ] )
        return results


# A function that calculates J-F's magic Tmunu observables (~ int dV T^{mu,nu} )
# Given a y-cut = [-1,1], and a pT-cut = [.2,5.0] e.g.
def calc_Tmunu(ET, y, pT, phi, ycut, pTcut):
        Tmunu = np.zeros(10)

        cut = (pTcut[0]<pT) & (pT<pTcut[1]) & (ycut[0]<y) & (y<ycut[1])
        # apply cut
        y = y[cut]
        ET = ET[cut]
        phi = phi[cut]
        pT = pT[cut]

        Pmu = np.array( [ET * np.cosh(y), pT * np.cos(phi), pT * np.sin(phi), ET * np.sinh(y) ] )

        # T0nu
        for i in range(4):
                Tmunu[i] = Pmu[i].sum()
        # Tij
        n = 4
        for i in range(1,4):
                for j in range(i, 4):
                        Tmunu[n] = (Pmu[i]*Pmu[j]/Pmu[0]).sum()
                        n += 1

        return Tmunu

# loop over the datatype iteratively and print every thing for a check
def dprint(data, dtype, prefix):
        for (name, sub_dtype, _) in dtype:
                sub_prefix = prefix+"/"+name
                if type(sub_dtype) == list:
                        dprint(data[name], sub_dtype, sub_prefix)
                else:
                        print(sub_prefix,'=', data[name][0])


def file_len(fname):
        """ do line count """
        p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result, err = p.communicate()
        if p.returncode != 0:
                raise IOError(err)
        return int(result.strip().split()[0])

def run_cmd(*args):
        """
        Run and log a subprocess.

        """
        cmd = ' '.join(args)
        logging.info('running command: %s', cmd)

        try:
                proc = subprocess.run(
                        cmd.split(), check=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        #stdout=sys.stdout, stderr=sys.stderr,
                        universal_newlines=True
                )
        except subprocess.CalledProcessError as e:
                logging.error(
                        'command failed with status %d:\n%s',
                        e.returncode, e.output.strip('\n')
                )
                raise
        else:
                logging.debug(
                        'command completed successfully:\n%s',
                        proc.stdout.strip('\n')
                )
                return proc


class Parser(argparse.ArgumentParser):
        """
        ArgumentParser that parses files with 'key = value' lines.

        """
        def __init__(self, *args, fromfile_prefix_chars='@', **kwargs):
                super().__init__(
                        *args, fromfile_prefix_chars=fromfile_prefix_chars, **kwargs
                )

        def convert_arg_line_to_args(self, arg_line):
                # split each line on = and prepend prefix chars to first arg so it is
                # parsed as a long option
                args = [i.strip() for i in arg_line.split('=', maxsplit=1)]
                args[0] = 2*self.prefix_chars[0] + args[0]
                return args


parser = Parser(
        usage=''.join('\n  %(prog)s ' + i for i in [
                '[options] <results_file>',
                'checkpoint <checkpoint_file>',
                '-h | --help',
        ]),
        description='''
Run relativistic heavy-ion collision events.

In the first form, run events according to the given options (below) and write
results to binary file <results_file>.

In the second form, run the event saved in <checkpoint_file>, previously
created by using the --checkpoint option and interrupting an event in progress.
''',
        formatter_class=argparse.RawDescriptionHelpFormatter
)


def parse_args():
        """
        Parse command line arguments according to the parser usage info.  Return a
        tuple (args, ic) where `args` is a normal argparse.Namespace and `ic` is
        either None or an np.array of the checkpointed initial condition.

        First, check for the special "checkpoint" form, and if found, load and
        return the args and checkpoint initial condition from the specified file.
        If not, let the parser object handle everything.

        This is a little hacky but it works fine.  Truth is, argparse can't do
        exactly what I want here.  I suppose `docopt` might be a better option, but
        it's not worth the effort to rewrite everything.

        """
        def usage():
                parser.print_usage(sys.stderr)
                sys.exit(2)

        if len(sys.argv) == 1:
                usage()

        return parser.parse_args()


parser.add_argument(
        'results', type=os.path.abspath,
        help=argparse.SUPPRESS
)
parser.add_argument(
        '--nevents', type=int, metavar='INT',
        help='number of events to run (default: run until interrupted)'
)
#parser.add_argument(
#       '--task', type=int, metavar='INT',
#       help='pass task id here'
#)
parser.add_argument(
        '--rankvar', metavar='VAR',
        help='environment variable containing process rank'
)
parser.add_argument(
        '--rankfmt', metavar='FMT',
        help='format string for rank integer'
)
parser.add_argument(
        '--tmpdir', type=os.path.abspath, metavar='PATH',
        help='temporary directory (default: {})'.format(tempfile.gettempdir())
)
parser.add_argument(
        '--startdir', type=os.path.abspath, metavar='PATH',
        help='directory where input parameter files are located'
)
parser.add_argument(
        '--tablesdir', type=os.path.abspath, metavar='PATH',
        help='directory where tables for iS3D, smash, music are located'
)
parser.add_argument(
        '--logfile', type=os.path.abspath, metavar='PATH',
        help='log file (default: stdout)'
)
parser.add_argument(
        '--loglevel', choices={'debug', 'info', 'warning', 'error', 'critical'},
        default='info',
        help='log level (default: %(default)s)'
)

class StopEvent(Exception):
        """ Raise to end an event early. """


#useful to keep track of the maximum radius of freezeout surface
def find_fireball_edge(surface_file):
        #open the freezeout surface file (MUSIC format)
        # tau, x, y, eta, ds0, ds1, ds2, ds3, ....
        #surf_tau = pd.read_csv( 'surface.dat', header=None, delimiter=' ', usecols=[0] )
        surf_x = pd.read_csv( surface_file, header=None, delimiter=' ', usecols=[1] , names=['surf_x'])
        surf_y = pd.read_csv( surface_file, header=None, delimiter=' ', usecols=[2] , names=['surf_y'])
        #surf_eta = pd.read_csv( 'surface.dat', header=None, delimiter=' ', usecols=[3] )
        surf_x_list = surf_x['surf_x'].tolist()
        surf_y_list = surf_y['surf_y'].tolist()
        x_max = 0.0
        y_max = 0.0
        for i in range( 0, len(surf_x_list) ):
                x_max = max( x_max, abs( surf_x_list[i]) )
        for i in range( 0, len(surf_y_list) ):
                y_max = max( y_max, abs( surf_y_list[i]) )
        x_y_max = max(x_max, y_max)
        return x_y_max

#where the magic happens
def run_events(args, results_file):
        """
        Run events as determined by user input:

                - Read options from `args`, as returned by `parser.parse_args()`.
                - Write results to binary file object `results_file`.

        Return True if at least one event completed successfully, otherwise False.

        """
        results = np.empty((), dtype=result_dtype)

        Nevents = args.nevents
        def run_single_event(event_number):
                """ run a single jetscape event"""
                results.fill(0)
                start_h = datetime.datetime.now()
                logging.info('TRENTO_FS_HYDRO started at %s ', start_h)
                #os.mkdir('output') #make a directory for freestream-milne to dump initial T00
                h = run_cmd('./TRENTO_FS_HYDRO') # this executable runs TRENTo, freestream-milne and MUSIC. writes fo surf to surface.dat
                #eps = calc_ecc() #call calc_ecc to calculate eccentricities of initial condition
                #print(h.stdout)
                #print(h.stderr)
                end_h = datetime.datetime.now()
                logging.info('TRENTO_FS_HYDRO finished in %s ', end_h - start_h)
                #copy the eccentricity files to final destination
                #subprocess.call( "mkdir {:s}/eccentricities".format( os.path.dirname(args.results) ) , shell=True)
                #for ecc_file in ['ecc_before_fs', 'ecc_r3_before_fs', 'ecc_r5_before_fs', 'ecc_before_hydro', 'ecc_r3_before_hydro', 'ecc_r5_before_hydro']:
                #        subprocess.call("mv {:s}.dat {:s}/eccentricities/{:s}-{:s}".format(ecc_file, os.path.dirname(args.results),
                #                                                                        ecc_file, os.path.basename(args.results), ), shell=True)
                # these files contain instances of viscous regulation inside high temperature region
                for rf, short in zip(['music_inv_reynolds_shear.txt', 'music_inv_reynolds_bulk.txt'], ['iR-shear','iR-bulk']):
                        if os.path.exists(rf):
                                logging.info("{:s} exists".format(rf))
                                subprocess.call("mv {:s} {:s}/{:s}-{:s}".format(rf, os.path.dirname(args.results), short, os.path.basename(args.results), ), shell=True)
                if os.stat("surface.dat").st_size == 0:
                        raise StopEvent('empty surface')
                else:
                        subprocess.call("mv surface.dat input/surface.dat", shell=True)

                #logging.info('Maximum radius of freezeout surface : %f fm', find_fireball_edge('input/surface.dat') )
                logging.info('Max radius of freezeout surface : {:f}'.format( find_fireball_edge('input/surface.dat') ) )
                # loop over 'delta-f'. in iS3D, delta_f = 1 : 14 mom, 2 : RTA C.E., 3 : McNelis Mod, 4 : Bernhard Mod
                # idf = 0,1,2,3 runs iS3D and smash. idf = 4 runs iS3D and urqmd-afterburner
                #for idf in range(0,4):
                for idf in [0]:
                        logging.info('idf  = {:d}'.format(idf))

                        if (idf < 4):
                                #use appropriate delta-f mode and smash box particle list
                                subprocess.call("sed -i \"s/df_mode.*/df_mode = {:d}/g\" iS3D_parameters.dat".format(idf+1), shell=True)
                                subprocess.call("sed -i \"s/hrg_eos.*/hrg_eos = {:d}/g\" iS3D_parameters.dat".format(3), shell=True)
                                subprocess.call("cp PDG/chosen_particles_box.dat PDG/chosen_particles.dat", shell=True)

                                subprocess.call("rm -f particles.dat", shell=True)
                                nparts = 0
                                logging.info("Run iS3D and smash")
                                start_s = datetime.datetime.now()
                                logging.info('SAMPLER_AFTERBURNER started at %s ', start_s)
                                sample_burn = run_cmd("./SAMPLER_AFTERBURNER") # this executable reads surface.dat into iS3D, samples, runs smash
                                #print(sample_burn.stdout)
                                #print(sample_burn.stderr)
                                end_s = datetime.datetime.now()
                                logging.info('SAMPLER_AFTERBURNER finished in %s ', end_s - start_s)
                                #split the file into iS3D hadrons and smash hadrons separately)
                                subprocess.call("csplit --digits=2  --quiet --prefix=hadron_list final_smash_hadrons.dat \"/# JetScape module: SMASH/+1\" \"{*}\" ", shell=True)
                                nparts += file_len("./hadron_list01") #particles after smash finished
                                subprocess.call("cat ./hadron_list01 >> particles.dat", shell=True)
                                #need to get the number of oversamples from first header line in final_smash_hadrons.dat
                                nos = 0
                                with open('final_smash_hadrons.dat') as f:
                                        first_line = f.readline().split()
                                        nos = int(first_line[1])
                                logging.info("{:d} particles in {:d} oversamples".format(nparts, nos))
                                # read final particle data after smash finished
                                with open("particles.dat", 'r') as f:
                                        particles = np.fromiter( (tuple(l.split()[2:]) for l in f if not l.startswith('#')), dtype=particle_dtype )

                        else:
                                #use delta_f mode 4 (Bernhard) and urqmd particle list
                                subprocess.call("sed -i \"s/df_mode.*/df_mode = {:d}/g\" iS3D_parameters.dat".format(4), shell=True)
                                subprocess.call("sed -i \"s/hrg_eos.*/hrg_eos = {:d}/g\" iS3D_parameters.dat".format(1), shell=True)
                                subprocess.call("cp PDG/chosen_particles_urqmd_v3.3+.dat PDG/chosen_particles.dat", shell=True)

                                subprocess.call("rm -f particles.dat", shell=True)
                                nparts = 0
                                logging.info("Run iS3D and urqmd-afterburner")
                                subprocess.call("mkdir results", shell=True)
                                start_s = datetime.datetime.now()
                                sample = run_cmd("iS3DSamplingOnly")
                                #print(sample.stdout)
                                #print(sample.stderr)
                                burn = run_cmd("afterburner results/particle_list_osc.dat particles.dat")
                                #print(burn.stdout)
                                #print(burn.stderr)
                                end_s = datetime.datetime.now()
                                logging.info('iS3DSamplingOnly and afterburner finished in %s ', end_s - start_s)
                                subprocess.call("rm -r -f results", shell=True)
                                nparts += file_len("./particles.dat")
                                #need to get the number of oversamples from first header line in final_iS3D_hadrons.dat
                                nos = 0
                                with open('final_iS3D_hadrons.dat') as f:
                                        first_line = f.readline().split()
                                        nos = int(first_line[1])
                                logging.info("{:d} particles in {:d} oversamples".format(nparts, nos))
                                #read particle data after urqmd finished
                                with open("particles.dat", 'r') as f:
                                        particles = np.fromiter( (tuple(l.split()) for l in f if not l.startswith('#')), dtype=particle_dtype )

                        logging.info('All Oversampling and Afterburners finished in %s ', end_s - end_h)
                        time_hydro = end_h - start_h
                        time_afterburn = end_s - end_h
                        total_time = end_s - start_h
                        logging.info('Fraction of time spent in trento, fs, hydro       : %f ', time_hydro / total_time)
                        logging.info('Fraction of time spent in sampling, afterburning  : %f ', time_afterburn / total_time)
                        logging.info('computing observables for delta_f = {:d}'.format(idf+1))

                        start_o = datetime.datetime.now()

                        # Save number of oversamples
                        results[expt_type][idf]['nsamples'] = nos

                        # Get pT, phi, eta, y ...
                        ET = particles['ET']
                        pT = particles['pT']
                        phi = particles['phi']
                        eta = particles['eta']
                        y = particles['y']
                        charge = particles['charge']
                        abs_pid = np.abs(particles['ID'])
                        abs_eta = np.fabs(eta)
                        abs_y = np.fabs(y)
                        
                        #TODO 
                        #NEED TO IMPLEMENT CUTS FOR DIFFERENT EXPERIMENTS/SYSTEMS
 
                        #STAR CUTS
                        #according to http://arxiv.org/pdf/0808.2041v2.pdf, for mean pT and dN/dy , pT is extrapolated using blast wave fit
                        #so no cuts for additional cuts for dN/dy or mean pT
                        #according to http://arxiv.org/pdf/1301.2187.pdf and http://arxiv.org/pdf/nucl-ex/0409033.pdf,
                        # v_n{k} are calculated with cuts : 0.15 < pT < 2.0 GeV/c ; |eta| < 1.0 
                        pT_cut_STAR  = [ 0.15, 2.0 ]
                        eta_cut_STAR = [-1.0, 1.0 ]
                        #STAR_kinematic_cuts = (eta > eta_cut_STAR[0]) & (eta < eta_cut_STAR[1]) & (pT > pT_cut_star[0]) & (pT < pT_cut_star[1])

                        # charged particle cut
                        charged = (charge != 0)
                        
                        #WE DONT HAVE THESE DATA FOR STAR
                        # 1) calculate dNch/eta and dET/deta
                        results[expt_type][idf]['dNch_deta'] = ( charged & (abs_eta < .5) ).sum() / (2 * .5) / nos
                        logging.info( "charged and abs_eta < 0.5 : " + str( (charged & (abs_eta < .5)).sum() / (2*.5) ) )
                        results[expt_type][idf]['dET_deta'] = ET[abs_eta < .6].sum() / (2 * .6) / nos
                        
                        #WE HAVE THESE DATA FOR STAR
                        # 2) calculate identified particle yield and mean_pT
                        for name, pid in species:
                                cut = (abs_pid == pid) & (abs_y < 0.5)
                                N = cut.sum()
                                results[expt_type][idf]['dN_dy'][name] = N * 1. / nos
                                results[expt_type][idf]['mean_pT'][name] = 0 if N==0 else pT[cut].mean()

                        #WE DONT HAVE THESE DATA FOR STAR
                        # 3) calculate chagred particle pT fluctuation
                        x = pT[charged & (abs_eta < .8) & (.15 < pT) & (pT < 2.)]
                        N = x.size
                        results[expt_type][idf]['pT_fluct_chg']['N'] = N
                        results[expt_type][idf]['pT_fluct_chg']['sum_pT'] = 0. if N==0 else x.sum()
                        results[expt_type][idf]['pT_fluct_chg']['sum_pT2'] = 0. if N==0 else (x**2).sum()

                        #WE DONT HAVE THESE DATA FOR STAR
                        # 4) calculate identified particle pT fluctuation
                        for name, pid in pi_K_p:
                                x = pT[(abs_pid == pid) & (abs_y < 0.5)]
                                N = x.size
                                results[expt_type][idf]['pT_fluct_pid'][name]['N'] = N
                                results[expt_type][idf]['pT_fluct_pid'][name]['sum_pT'] = 0. if N == 0 else x.sum()
                                results[expt_type][idf]['pT_fluct_pid'][name]['sum_pT2'] = 0. if N == 0 else (x**2).sum()
                        
                        #WE HAVE THESE DATA FOR STAR
                        # 5) pT-integrated Qvector with STAR cut
                        pTbins = np.array([ pT_cut_STAR ])
                        etabins = np.array([ eta_cut_STAR ])
                        info = calc_Qvector(pT[charged], phi[charged], eta[charged], pTbins, etabins, Nharmonic)
                        results[expt_type][idf]['flow']['N'] = info['N'][0,0]
                        results[expt_type][idf]['flow']['Qn'] = info['Qn'][0,0,:]

                        #WE DONT HAVE THESE DATA FOR STAR
                        # 6) The magic Tmunu
                        results[expt_type][idf]['Tmunu'] = calc_Tmunu(ET, y, pT, phi, ycut=[-0.6, 0.6], pTcut=[0., 10.]) / nos
                        results[expt_type][idf]['Tmunu_chg'] = calc_Tmunu(ET[charged], y[charged], pT[charged], phi[charged], ycut=[-0.6, 0.6], pTcut=[0., 10.]) / nos

                        #WE DONT HAVE THESE DATA FOR STAR
                        # 7) Identified differential flow for charged stable hadrons
                        pTbins=np.array([[Qn_diff_pT_cuts[i],Qn_diff_pT_cuts[i+1]] for i in range(Qn_diff_NpT)])
                        rapbins = np.array([[-1, 1]])
                        for name, pid in Qn_species:
                                if name == 'Ch':
                                    cut = charged
                                else:
                                    cut = (abs_pid == pid)
                                info = calc_Qvector(pT[cut], phi[cut], y[cut], pTbins, rapbins, Nharmonic_diff)
                                results['d_flow_pid'][idf][name]['N'] = info['N'][:,0]
                                results['d_flow_pid'][idf][name]['Qn'] = info['Qn'][:,0,:]

                        end_o = datetime.datetime.now()
                        logging.info('Reading hadron list and computing observables finished in %s ', end_o - start_o)

        nfail = 0
        # run each initial condition event and save results to file
        for n in range(Nevents):
                logging.info('starting jetscape event %d', n)
                try:
                        run_single_event(n)
                except StopEvent as e:
                        logging.info('event stopped: %s', e)
                except Exception:
                        logging.exception('event %d failed', n)
                        nfail += 1
                        if nfail/Nevents > .25:
                                logging.critical('too many failures, stopping events')
                                break
                        logging.warning('continuing to next event')
                        continue

                results_file.write(results.tobytes())
                logging.info('event %d completed successfully', n)

        return nfail


def main():
        args = parse_args()

        # starting fresh -> truncate output files
        filemode = 'w'

        # must handle rank first since it affects paths
        if args.rankvar:
                rank = os.getenv(args.rankvar)
                if rank is None:
                        sys.exit('rank variable {} is not set'.format(args.rankvar))

                if args.rankfmt:
                        rank = args.rankfmt.format(int(rank))

                # append rank to path arguments, e.g.:
                #   /path/to/output.log  ->  /path/to/output/<rank>.log
                for a in ['results', 'logfile']:
                        value = getattr(args, a)
                        if value is not None:
                                root, ext = os.path.splitext(value)
                                setattr(args, a, os.path.join(root, rank) + ext)


        os.makedirs(os.path.dirname(args.results), exist_ok=True)

        if args.logfile is None:
                logfile_kwargs = dict(stream=sys.stdout)
        else:
                logfile_kwargs = dict(filename=args.logfile, filemode=filemode)
                os.makedirs(os.path.dirname(args.logfile), exist_ok=True)

        logging.basicConfig(
                level=getattr(logging, args.loglevel.upper()),
                format='[%(levelname)s@%(relativeCreated)d] %(message)s',
                **logfile_kwargs
        )
        logging.captureWarnings(True)

        start = datetime.datetime.now()

        logging.info('started at %s', start)
        logging.info('arguments: %r', args)

        # translate SIGTERM to KeyboardInterrupt
        signal.signal(signal.SIGTERM, signal.default_int_handler)
        logging.debug('set SIGTERM handler')
        #prefix=str(args.startdir)
        tablesdir=str(args.tablesdir)
        logging.info('starting dir: {:s}'.format(tablesdir))
        with open(args.results, filemode + 'b') as results_file, tempfile.TemporaryDirectory( prefix='jetscape-', dir=args.tmpdir) as workdir:
                os.chdir(workdir)
                os.mkdir("input/")
                #note in general symlinks are not safe if these files need to be modified!
                #copy input files from corresponding design point directory
                os.system( 'cp {:s}/jetscape_init.xml jetscape_init.xml'.format(tablesdir)  )
                os.system( 'cp {:s}/freestream_input freestream_input'.format(tablesdir)  )
                os.system( 'cp {:s}/music_input music_input'.format(tablesdir) )
                os.system( 'cp {:s}/iS3D_parameters.dat iS3D_parameters.dat'.format(tablesdir)  )
                #get common tables necessary for iS3D, smash, music
                os.system( 'cp -R {:s}/PDG PDG'.format(tablesdir)  )
                os.system( 'cp -R {:s}/EOS EOS'.format(tablesdir)  )
                os.system( 'cp -R {:s}/tables tables'.format(tablesdir)  )
                os.system( 'cp -R {:s}/deltaf_coefficients deltaf_coefficients'.format(tablesdir)  )
                os.system( 'cp -R {:s}/smash_input smash_input'.format(tablesdir) )
                #copy the executables to cwd
                os.system( 'cp {:s}/TRENTO_FS_HYDRO TRENTO_FS_HYDRO'.format(tablesdir)  )
                os.system( 'cp {:s}/SAMPLER_AFTERBURNER SAMPLER_AFTERBURNER'.format(tablesdir)  )
                

                logging.info('working directory: %s', workdir)

                try:
                        status = run_events(args, results_file) #run_events returns the number of failed events
                except KeyboardInterrupt:
                        # after catching the initial SIGTERM or interrupt, ignore them
                        # during shutdown -- this ensures everything will exit gracefully
                        # in case of additional signals (short of SIGKILL)
                        signal.signal(signal.SIGTERM, signal.SIG_IGN)
                        signal.signal(signal.SIGINT, signal.SIG_IGN)
                        status = True
                        logging.info( 'interrupt or signal at %s, cleaning up...', datetime.datetime.now() )
                        if args.checkpoint is not None:
                                logging.info( 'current event saved in checkpoint file %s', args.checkpoint )

        end = datetime.datetime.now()
        logging.info('finished at %s, %s elapsed', end, end - start)

        #if events fail or keyboard interrupt status > 0
        if status:
                sys.exit(1)


if __name__ == "__main__":
        main()
