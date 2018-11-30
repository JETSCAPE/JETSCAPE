#!/usr/bin/env python3

import argparse
import datetime
from itertools import chain, repeat
import logging
import math
import os
import pickle
import signal
import subprocess
import sys
import tempfile
import numpy as np

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
]
# fully specify numeric data types, including endianness and size, to
# ensure consistency across all machines
float_t = '<f8'
int_t = '<i8'
complex_t = '<c16'
# results "array" (one element)
# to be overwritten for each event
result_data_type=[
    ('initial_entropy', float_t), # 1) initial entropy (eta=0) of TRENTo
    ('dNch_deta', float_t),       # 2) dNch/deta at eta=0
    ('dET_deta', float_t),        # 3) dET/deta at eta=0
    ('dN_dy', [(s, float_t) 
        for (s, _) in species]),  # 4) identified particle yield
    ('mean_pT', [(s, float_t) 
        for (s, _) in species]),  # 5) identified particle <pT> 
    ('mean_pT2', [(s, float_t) 
        for (s, _) in species]),  # 6) identified particle <pT^2> 
    ('pT_fluct', [('N', int_t), ('sum_pT', float_t), ('sum_pTsq', float_t)]),
                                  # 7) pT fluct = sum_pTsq/N^2-(sum_pT/N)^2
    ('flow', [('N', int_t), ('Qn', complex_t, 8)]),
                                  # 8) pT-integrated Q-vector
                                  #    of charged particles,
    ('EP', [('Vn', float_t, 8), ('Psin', float_t, 8)]),
                                  # 9) Event-plane of charged particles,
                                  #    for comparison purpose
    ('d-flow', [(s, [('N', int_t, [6, 3]), ('Qn', complex_t, [6, 3, 4])])
        for s in ['pion','kaon','proton']]), 
                                  # 10) Double differential q-vector
                                  #    of proton, kaon and pions
                                  #    A grid of 10 pT-points * 7 eta-points
    ('d-ref', [('N', int_t), ('Qn', complex_t, 4)]),
                                  # 11) Q-vector for 
                                  #    of charged particles,
]


def calculate_obs(input_filename, output_filename):
    # A function that calculates event-planes given a list of azimuthal angles
    def calculate_eventplanes(phi):
        Vn_table = {n: float_t for n in range(1,8)}
        Psin_table =  {n: float_t for n in range(1,8)}
        N = phi.size
        for n in range(1,8):
            qy = np.sin(n*phi) # array of cos(n*phi_i)
            Qy = qy.sum()       # Qx = N<cos(n*phi)>
            qx = np.cos(n*phi) # array of cos(n*phi_i)
            Qx = qx.sum()       # Qy = N<sin(n*phi)>
            cov = np.cov([qy*N, qx*N])/(N-1.) # covariance between Qy and Qx array
            varQy = cov[0,0]; varQx = cov[1,1]; covQxQy = cov[0,1]
            Psin = np.arctan2(Qy, Qx)/n
            Psinerr = 1./n*np.sqrt(varQy*Qx**2 + varQx*Qy**2 - 2.*covQxQy*Qy*Qx) \
                  /(Qy**2+Qx**2) # uncertainty of Psin
            Vn = np.sqrt((Qx**2+Qy**2 - N)/N/(N-1.))
            Vnerr = np.sqrt(Qx**2*varQx + Qy**2*varQy + 2.*Qx*Qy*covQxQy)\
                        /N/(N - 1.)/Vn # uncertainty of Vn
            Vn_table[n] = Vn
            Psin_table[n] = Psin
        return Vn_table, Psin_table
   
    # A function that calculates double differential q-vectos
    # given a list of pT-bins and a list of etabins
    # for example, pTbins = np.array([[0,1],[1,2],[2,3]])
    # and etabins = np.array([[-2.5,-1.5],[-1.5,-.5],[-.5,.5],[.5,1.5],[1.5,2.5]])
    def calc_Qvector(pT, phi, eta, pTbins, etabins, order=8):
        grid_shape = [pTbins.shape[0], etabins.shape[0]]
        # calculate Qn
        results = {'N': np.zeros(grid_shape, dtype=int_t),
                   'Qn': np.zeros([*grid_shape, order], dtype=complex_t)}
        for i, (pT_low, pT_high) in enumerate(pTbins):
            pT_cut = (pT_low<=pT) & (pT<=pT_high) 
            for j, (eta_low, eta_high) in enumerate(etabins):
                eta_cut = (eta_low<=eta) & (eta<=eta_high)        
    
                kinematic_cut = pT_cut & eta_cut
                # apply cut
                dphi = phi[kinematic_cut]
                results['N'][i,j] = dphi.size
                results['Qn'][i,j] = np.array([ np.sum(np.exp(1j*n*dphi)) 
                                                for n in range(1, 1+order) ])
        return results

    # read final particle data (output from afterburner)
    with open(input_filename, 'r') as f:
        parts = np.fromiter(
            (tuple(l.split()) for l in f if not l.startswith('#')),
            dtype=[
                ('ID', int),
                ('charge', int),
                ('pT', float),
                ('ET', float),
                ('mT', float),
                ('phi', float),
                ('y', float),
                ('eta', float),
            ]
        )

    logging.info('computing observables')
    results = np.empty((), dtype=result_data_type)
    results.fill(0.) # initialize all results entry with 0

    pT = parts['pT']
    phi = parts['phi']
    eta = parts['eta']
    y = parts['y']
    charge = parts['charge']
    abs_ID = np.abs(parts['ID'])
    abs_eta = np.fabs(eta)

    # some useful cut
    charged = (charge != 0)
    midrapidity = (np.fabs(y) < .5)

    # calcualte dNch/eta and dET/deta
    results['dNch_deta'] = np.count_nonzero(charged & (abs_eta < .5)) / 1.
    results['dET_deta'] = parts['ET'][abs_eta < .6].sum() / (2*.6)

    # identified particle yield
    for name, i in species:
        cut = (abs_ID == i) & midrapidity
        N = np.count_nonzero(cut)
        results['dN_dy'][name] = N
        results['mean_pT'][name] = (0. if N == 0 else pT[cut].mean())
        pTsq = pT[cut] * pT[cut]
        results['mean_pT2'][name] = (0. if N == 0 else pTsq.mean())

    # pT fluctuation
    pT_alice = pT[charged & (abs_eta < .8) & (.15 < pT) & (pT < 2.)]
    results['pT_fluct']['N'] = pT_alice.size
    results['pT_fluct']['sum_pT'] = pT_alice.sum()
    results['pT_fluct']['sum_pTsq'] = np.inner(pT_alice, pT_alice)

    # pT-integrated Qvector with ALICE cut 
    pTbins = np.array([[.2, 5.0]])
    etabins = np.array([[-1,1]])
    Q = calc_Qvector(pT[charged], phi[charged], eta[charged], 
                     pTbins, etabins, order=8)
    results['flow']['N'] = Q['N'][0,0]

    results['flow']['Qn'] = Q['Qn'][0,0,:]

    # identified particle pT-differential Qvector and reference Qvector
    pTbins = np.array([[0,.5],[.5,1.],[1.,1.5],[1.5,2.],[2.,2.5],[2.5,3.]])
    etabins = np.array([[-3,-1],[-1,1],[1, 3]])
    ref_pTbins = np.array([[0.2,5.]])
    ref_etabins = np.array([[-2,2]])
    for s, pid in zip(['pion','kaon','proton'], [211, 321, 2212]):
        cut = (abs_ID == pid) & charged
        Q = calc_Qvector(pT[cut], phi[cut], eta[cut], 
                         pTbins, etabins, order=4)
        Qref = calc_Qvector(pT[charged], phi[charged], eta[charged], 
                         ref_pTbins, ref_etabins, order=4)
        results['d-flow'][s]['N'] = Q['N']
        results['d-flow'][s]['Qn'] = Q['Qn']
        results['d-ref']['N'] = Qref['N'][0,0]
        results['d-ref']['Qn'] = Qref['Qn'][0,0,:]

    # write results to file
    with open(output_filename, 'wb') as results_file:
        results_file.write(results.tobytes())
   
def view_results(filename):
    data = np.fromfile(filename, dtype=result_data_type)
    print("pT-integrate Qvec:N:", *data['flow']['N'])
    print("pT-integrate Qvec:Qn:", *data['flow']['Qn'])
    print("double-diff Qvec:pion:N:", *data['d-flow']['pion']['N'])
    print("double-diff Qvec:pion:Qn:", *data['d-flow']['pion']['Qn'])

def main():
    calculate_obs("final_smash_hadrons.dat", "results.dat")
    view_results("results.dat")

if __name__ == "__main__":
    main()
