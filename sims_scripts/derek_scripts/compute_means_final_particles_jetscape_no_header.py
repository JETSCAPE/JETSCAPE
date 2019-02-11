import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import math
import sys
#from statistics import *

#particle_list = pd.read_csv( sys.argv[1], sep=' ', names=['idx', 'H','status', 'mcid', 'label', 'pT', 'eta', 'phi', 'e', 't', 'x', 'y', 'z'] )

#species dependent info
#pion 211
pi_pT  = []
pi_eta = []
pi_phi = []
pi_y = []

#kaon 321
k_pT  = []
k_eta = []
k_phi = []
k_y = []

#proton 2212
p_pT  = []
p_eta = []
p_phi = []
p_y = []

#all particles
all_t = []
all_x = []

#after cut on pseudorap or rap
pi_pT_mid = []
k_pT_mid = []
p_pT_mid = []

pi_phi_mid = []
k_phi_mid = []
p_phi_mid = []

all_t_mid = []
all_x_mid = []

sample_dir = sys.argv[1]
nsamples = int(sys.argv[2])

#the range of psedoraprapidity cut
etamax = 1.0

#the range of rapidity cut
ymax = 1.0

for sample in range(0, nsamples):
    hadron_file = sample_dir + 'sample_' + str(sample) + '/finaliSSHadrons.dat'
    print("Reading file : " + hadron_file )
    mcid_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[3], skiprows=6, names=['mcid_df'])
    mcid = mcid_df['mcid_df'].tolist()
    pT_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[5], skiprows=6, names=['pT_df'])
    pT = pT_df['pT_df'].tolist()
    eta_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[6], skiprows=6, names=['eta_df'])
    eta = eta_df['eta_df'].tolist()
    phi_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[7], skiprows=6, names=['phi_df'])
    phi = phi_df['phi_df'].tolist()
    e_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[8], skiprows=6, names=['e_df'])
    e = e_df['e_df'].tolist()
    t_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[9], skiprows=6, names=['t_df'])
    t = t_df['t_df'].tolist()
    x_df = pd.read_csv( hadron_file, header=None, delimiter=' ', usecols=[10], skiprows=6, names=['x_df'])
    x = x_df['x_df'].tolist()

    #3d momentum space lists
    for i in range(0, len(e) ):
        if ( mcid[i] == 211 ):
            pi_pT.append( pT[i] )
            pi_eta.append( eta[i] )
            pi_phi.append( phi[i] )
            pz = pT[i] * math.sinh( eta[i] )
            energy = e[i]
            arg = (energy + pz) / (energy - pz)
            if (pz == 0):
                y = 0
            else:
                y = 0.5 * math.log( abs(arg) )
            pi_y.append( y )

        elif ( mcid[i] == 321 ):
            k_pT.append( pT[i] )
            k_eta.append( eta[i] )
            k_phi.append( phi[i] )
            pz = pT[i] * math.sinh( eta[i] )
            energy = e[i]
            arg = (energy + pz) / (energy - pz)
            if (pz == 0):
                y = 0
            else:
                y = 0.5 * math.log( abs(arg) )
            k_y.append( y )

        elif ( mcid[i] == 2212 ):
            p_pT.append( pT[i] )
            p_eta.append( eta[i] )
            p_phi.append( phi[i] )
            pz = pT[i] * math.sinh( eta[i] )
            energy = e[i]
            arg = (energy + pz) / (energy - pz)
            if (pz == 0):
                y = 0
            else:
                y = 0.5 * math.log( abs(arg) )
            p_y.append( y )
    #spacetime info with pseudorap cut
    for i in range(0, len(t)):
        pz = pT[i] * math.sinh( eta[i] )
        energy = e[i]
        arg = (energy + pz) / (energy - pz)
        if (pz == 0):
            y = 0
        else:
            y = 0.5 * math.log( abs(arg) )
        if ( abs( y ) < ymax ):
            all_t_mid.append( t[i] )
            all_x_mid.append( x[i] )

#pseudorap cuts on momentum space info 
for i in range(0, len(pi_pT) ):
    if ( abs(pi_y[i]) < ymax ):
        pi_pT_mid.append( pi_pT[i] )
        pi_phi_mid.append( pi_phi[i] )
for i in range(0, len(k_pT) ):
    if ( abs(k_y[i]) < ymax ):
        k_pT_mid.append( k_pT[i] )
        k_phi_mid.append( k_phi[i] )
for i in range(0, len(p_pT) ):
    if ( abs(p_y[i]) < ymax ):
        p_pT_mid.append( p_pT[i] )
        p_phi_mid.append( p_phi[i] )

#compute mean values

#mean yield
mean_yield_pi_midrap = len(pi_pT)
mean_yield_k_midrap = len(k_pT)
mean_yield_p_midrap = len(p_pT)

#mean pT
mean_pT_pi_midrap = np.mean(pi_pT)
mean_pT_k_midrap = np.mean(k_pT)
mean_pT_p_midrap = np.mean(p_pT)

#pT width
var_pT_pi_midrap = np.var(pi_pT)
var_pT_k_midrap = np.var(k_pT)
var_pT_p_midrap = np.var(p_pT)

#write to files 
file = open("mean_vals.dat","w")
file.write("pi : ( N, <pT>, sigma_pT ) = ( " + str(mean_yield_pi_midrap) + ", " + str(mean_pT_pi_midrap) + ", " + str(var_pT_pi_midrap) + " ) \n" )
file.write("k  : ( N, <pT>, sigma_pT ) = ( " + str(mean_yield_k_midrap) + ", " +  str(mean_pT_k_midrap) + ", " + str(var_pT_k_midrap)  + " ) \n" )
file.write("p  : ( N, <pT>, sigma_pT ) = ( " + str(mean_yield_p_midrap) + ", " +  str(mean_pT_p_midrap) + ", " + str(var_pT_p_midrap)  + " ) \n" )
file.close()
