#!/usr/bin/env python3

import matplotlib.pyplot as plt
import logging
import math
import os
import sys
import numpy as np
logging.getLogger().setLevel(logging.DEBUG)

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
]
pi_K_p = [
	('pion', 211),
	('kaon', 321),
	('proton', 2212),
]

NpT = 10
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
result_dtype=[
('initial_entropy', float_t, 1), 
('impact_parameter', float_t, 1), 
('npart', float_t, 1), 
('ALICE',
	[
		# 1) dNch/deta, eta[-0.5, 0.5], charged
		('dNch_deta', float_t, 1),	   
		# 2) dET/deta, eta[-0.6, 0.6]
		('dET_deta', float_t, 1),		
		# 3.1) The Tmunu observables, eta[-0.6, 0.6]
		('Tmunu', float_t, 10),	 
		# 3.2) The Tmunu observables, eta[-0.5, 0.5], charged   
		('Tmunu_chg', float_t, 10),	
		# 4.1) identified particle yield
		('dN_dy', 	[(name, float_t, 1) for (name,_) in species], 1),  
		# 4.2) identified particle <pT>
		('mean_pT', [(name, float_t, 1) for (name,_) in species], 1),  
		# 5.1) pT fluct, pT[0.15, 2], eta[-0.8, 0.8], charged
		('pT_fluct_chg', [	('N', int_t, 1), 
							('sum_pT', float_t, 1), 
							('sum_pT2', float_t, 1)], 1), 
		# 5.2) pT fluct, pT[0.15, 2], eta[-0.8, 0.8], pi, K, p			
		('pT_fluct_pid', [	(name, [	('N', int_t, 1), 
										('sum_pT', float_t, 1), 
										('sum_pT2', float_t, 1)], 1	)
							  for (name,_) in pi_K_p	], 1), 
		# 6) Q vector, pT[0.2, 5.0], eta [-0.8, 0.8], charged
		('flow', [	('N', int_t, 1), 
					('Qn', complex_t, Nharmonic)], 1),
		# 7) Q vector, diff-flow eta[-0.8, 0.8], pi, K, p	
		# It uses #6 as its reference Q vector		
		('d_flow_chg', [('N', int_t, NpT), 
						('Qn', complex_t, [NpT, Nharmonic_diff])], 1), 
		('d_flow_pid', [(name, [('N', int_t, NpT), 
								('Qn', complex_t, [NpT, Nharmonic_diff])], 1)
						for (name,_) in pi_K_p	], 1), 
	], 1)
]


# A function that calculates event-planes given a list of azimuthal angles
def calculate_eventplanes(phi, Nharmonic):
	V_n = {n: float_t for n in range(1,Nharmonic+1)}
	Psi_n =  {n: float_t for n in range(1,Nharmonic+1)}
	N = phi.size
	Npairs = N*(N-1.)
	for n in range(1,Nharmoic+1):
		Qx = np.cos(n*phi).sum() 
		Qy = np.sin(n*phi).sum()
		V_n[n] = np.sqrt((Qx**2+Qy**2 - N)/pairs)
		Psi_n[n] = np.arctan2(Qy, Qx)/n
	return V_n, Psi_n

# A function that calculates double differential q-vectors
# given a list of pT-bins and a list of eta-bins
# for example:
#	 pTbins = np.array([[0,1],[1,2],[2,3]])
#	etabins = np.array([[-3,-1],[-1,1],[1,3]])
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

			results['Qn'][i,j] = np.array([ np.exp(1j*n*sub_phi).sum()
										   for n in range(1, 1+Nharmonic) ])
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

	Pmu = np.array([ET * np.cosh(y), pT * np.cos(phi),
					pT * np.sin(phi), ET * np.sinh(y) ])
	
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

def calculate_obs(input_filename, output_filename):
	# read final particle data (output from afterburner)
	with open(input_filename, 'r') as f:
		particles = np.fromiter(
			(tuple(l.split()[2:]) for l in f if not l.startswith('#')),
			dtype=particle_dtype
		)

	logging.info('computing observables')

	# Create result data entry and set all to 0
	results = np.empty((), dtype=result_dtype)
	results.fill(0.) 

	# Get pT, phi, eta, y ...
	ET = particles['ET']
	pT = particles['pT']
	phi = particles['phi']
	eta = particles['eta']
	y = particles['y']
	charge = particles['charge']
	abs_pid = np.abs(particles['ID'])
	abs_eta = np.fabs(eta)
	abs_y = np.fabs(eta)

        #plot dN/dphi 
	#plt.hist(phi, 20)
	#plt.show()

	# charged particle cut
	charged = (charge != 0)

	# 1) calculate dNch/eta and dET/deta
	results['ALICE']['dNch_deta'] = (charged & (abs_eta < .5)).sum() / (2*.5)
	results['ALICE']['dET_deta'] = ET[abs_eta < .6].sum() / (2*.6)

	# 2) calculate identified particle yield and mean_pT
	for name, pid in species:
		cut = (abs_pid == pid) & (abs_y < 0.5)
		N = cut.sum()
		results['ALICE']['dN_dy'][name] = N
		results['ALICE']['mean_pT'][name] = 0 if N==0 else pT[cut].mean()

	# 3) calculate chagred particle pT fluctuation
	x = pT[charged & (abs_eta < .8) & (.15 < pT) & (pT < 2.)]
	N = x.size
	results['ALICE']['pT_fluct_chg']['N'] = N
	results['ALICE']['pT_fluct_chg']['sum_pT'] = 0. if N==0 else x.sum()
	results['ALICE']['pT_fluct_chg']['sum_pT2'] = 0. if N==0 else (x**2).sum()

	# 4) calculate identified particle pT fluctuation
	for name, pid in pi_K_p:
		x = pT[(abs_pid == pid) & (abs_y < 0.5)]
		N = x.size
		results['ALICE']['pT_fluct_pid'][name]['N'] = N
		results['ALICE']['pT_fluct_pid'][name]['sum_pT'] = 0. \
														if N == 0 else x.sum()
		results['ALICE']['pT_fluct_pid'][name]['sum_pT2'] = 0. \
														if N == 0 else (x**2).sum()

	# 5) pT-integrated Qvector with ALICE cut
	pTbins = np.array([[.2, 5.0]])
	etabins = np.array([[-0.8, 0.8]])
	info = calc_Qvector(pT[charged], phi[charged], eta[charged],
					 pTbins, etabins, Nharmonic)
	results['ALICE']['flow']['N'] = info['N'][0,0]
	results['ALICE']['flow']['Qn'] = info['Qn'][0,0,:]

	# 6) The magic Tmunu
	results['ALICE']['Tmunu'] = calc_Tmunu(ET, y, pT, phi, 
								ycut=[-0.6, 0.6], pTcut=[0., 10.])

	results['ALICE']['Tmunu_chg'] = calc_Tmunu(ET[charged], y[charged], 
								pT[charged], phi[charged], 
								ycut=[-0.6, 0.6], pTcut=[0., 10.])

	# 7) Diff flow
	pTbins = np.array([[0,.2], [.2,.4], [.4,.6],[.6,.8],[.8,1.],
			[1.,1.2], [1.2,1.5], [1.5,2.], [2.,2.5], [2.5,3]])
	etabins = np.array([[-0.8, 0.8]])
	info = calc_Qvector(pT[charged], phi[charged], eta[charged],
					 pTbins, etabins, Nharmonic_diff)
	results['ALICE']['d_flow_chg']['N'] = info['N'][:,0]
	results['ALICE']['d_flow_chg']['Qn'] = info['Qn'][:,0,:]
	for name, pid in pi_K_p:
		cut = (abs_pid == pid)
		info = calc_Qvector(pT[cut], phi[cut], eta[cut],
					 pTbins, etabins, Nharmonic_diff)
		results['ALICE']['d_flow_pid'][name]['N'] = info['N'][:,0]
		results['ALICE']['d_flow_pid'][name]['Qn'] = info['Qn'][:,0,:]

	# write results to file
	with open(output_filename, 'wb') as results_file:
		results_file.write(results.tobytes())

# loop over the datatype iteratively and print every thing for a check
def dprint(data, dtype, prefix):
	for (name, sub_dtype, _) in dtype:
		sub_prefix = prefix+"/"+name
		if type(sub_dtype) == list:
			dprint(data[name], sub_dtype, sub_prefix)
		else:
			print(sub_prefix,'=', data[name][0])

if __name__ == "__main__":
	calculate_obs("hadron_list01", "results.dat")
	data = np.fromfile("results.dat", dtype=result_dtype)
	dprint(data, result_dtype, '')
