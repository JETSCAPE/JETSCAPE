#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys, os, glob
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

def list2array(func):
	def func_wrapper(x, w):
		try:
			x = np.array(x)
			w = np.array(w)
		except:
			raise ValueError("cannot interpret input as numpy array...")
		return func(x, w)
	return func_wrapper

def weighted_mean_std(x, w=None):
	if w is None:
		Neff = x.size
		mean = np.mean(x)
		std = np.std(x)/np.sqrt(Neff-1.+1e-9)
	else:
		Neff = np.sum(w)**2/np.sum(w**2)
		mean = np.average(x, weights=w)
		std = ( np.average((x-mean)**2, weights=w**2)/(Neff-1.+1e-9) ) **.5
	return mean, std

def calculate_dNdeta(ds, exp, cen):
	Ne = len(ds)
	cenM = np.mean(cen, axis=1)
	index = (cen/100.*Ne).astype(int)
	obs = np.zeros_like(cenM)
	obs_err = np.zeros_like(cenM)
	for i, (nl, nh) in enumerate(zip(index[:,0], index[:,1])):
		obs[i], obs_err[i] = weighted_mean_std(ds[exp]['dNch_deta'][nl:nh])
	return {'Name': 'dNch_deta', 'cenM': cenM, 'pTM' : None,
			'obs': obs, 'err': obs_err}


def calculate_dETdeta(ds, exp, cen):
	Ne = len(ds)
	cenM = np.mean(cen, axis=1)
	index = (cen/100.*Ne).astype(int)
	obs = np.zeros_like(cenM)
	obs_err = np.zeros_like(cenM)
	for i, (nl, nh) in enumerate(zip(index[:,0], index[:,1])):
		obs[i], obs_err[i] = weighted_mean_std(ds[exp]['dET_deta'][nl:nh])
	return {'Name': 'dNch_deta', 'cenM': cenM, 'pTM' : None,
			'obs': obs, 'err': obs_err}
	
def calculate_dNdy(ds, exp, cen):
	Ne = len(ds)
	cenM = np.mean(cen, axis=1)
	index = (cen/100.*Ne).astype(int)
	obs = {s: np.zeros_like(cenM) for (s, _) in species}
	obs_err = {s: np.zeros_like(cenM) for (s, _) in species}
	for (s, _) in species:
		for i, (nl, nh) in enumerate(zip(index[:,0], index[:,1])):
			obs[s][i], obs_err[s][i] = weighted_mean_std(ds[exp]['dN_dy'][s][nl:nh])
	return {'Name': 'dNch_deta', 'cenM': cenM, 'pTM' : None,
			'obs': obs, 'err': obs_err}	
		
def calculate_mean_pT(ds, exp, cen):
	Ne = len(ds)
	cenM = np.mean(cen, axis=1)
	index = (cen/100.*Ne).astype(int)
	obs = {s: np.zeros_like(cenM) for (s, _) in species}
	obs_err = {s: np.zeros_like(cenM) for (s, _) in species}
	for (s, _) in species:
		for i, (nl, nh) in enumerate(zip(index[:,0], index[:,1])):
			obs[s][i], obs_err[s][i] = weighted_mean_std(ds[exp]['mean_pT'][s][nl:nh])
	return {'Name': 'dNch_deta', 'cenM': cenM, 'pTM' : None,
			'obs': obs, 'err': obs_err}
			
def calculate_mean_pT_fluct(ds, exp, cen):
	Ne = len(ds)
	cenM = np.mean(cen, axis=1)
	index = (cen/100.*Ne).astype(int)
	obs = np.zeros_like(cenM) 
	obs_err = np.zeros_like(cenM) 
	for (s, _) in species:
		for i, (nl, nh) in enumerate(zip(index[:,0], index[:,1])):
			N = ds[exp]['pT_fluct_chg']['N'][nl:nh]
			sumpT = ds[exp]['pT_fluct_chg']['sum_pT'][nl:nh]
			sumpT2 = ds[exp]['pT_fluct_chg']['sum_pT2'][nl:nh]
			fluct = sumpT2/N-(sumpT/N)**2
			weight = N**2
			obs[i], obs_err[i] = weighted_mean_std(fluct, weight)
	return {'Name': 'dNch_deta', 'cenM': cenM, 'pTM' : None,
			'obs': obs, 'err': obs_err}


def calculate_vn(ds, exp, cen):
	@list2array
	def obs_and_err(qn, m):
		N = len(m)
		w = m*(m-1.)
		cn2 = (np.abs(qn)**2 - m)/w
		avg_cn2, var_avg_cn2 = weighted_mean_std(cn2, w)
		vn = np.sqrt(avg_cn2)
		vn_err = np.sqrt(var_avg_cn2/4./avg_cn2)
		return vn, vn_err
	Ne = len(ds)
	cenM = np.mean(cen, axis=1)
	index = (cen/100.*Ne).astype(int)

	obs = np.zeros([len(cenM), Nharmonic])
	obs_err = np.zeros([len(cenM), Nharmonic])
	for i, (nl, nh) in enumerate(zip(index[:,0], index[:,1])):
		M = ds[exp]['flow']['N'][nl:nh] + 1e-15
		for n in range(Nharmonic): 
			Q = ds[exp]['flow']['Qn'][nl:nh, n]
			obs[i,n], obs_err[i,n] = obs_and_err(Q, M)
	return {'Name': 'vn', 'cenM': cenM, 'pTM' : None,
			'obs': obs, 'err': obs_err}

def calculate_diff_vn(ds, exp, cenbins, pTbins, pid='chg'):

	Ne = len(ds)
	pTbins = np.array(pTbins)
	cenbins = np.array(cenbins)
	cenM = np.mean(cenbins, axis=1)
	pTM = np.mean(pTbins, axis=1)
	Cindex = (cenbins/100.*Ne).astype(int)
	
	if pid == 'chg':
		obs = 'd_flow_chg'
		data = ds[exp][obs]
	else:
		obs = 'd_flow_pid'
		data = ds[exp][obs][s]
	
	# need soft flow within the same centrality bin first
	# only needs Ncen x [v2, v3]
	vnref = calculate_vn(ds, exp, cenbins)
	
	# calculate hard vn
	vn = np.zeros([len(cenM), len(pTM), Nharmonic_diff])
	vn_err = np.zeros([len(cenM), len(pTM), Nharmonic_diff])
	for i, (nl, nh) in enumerate(Cindex):
		for j, (pl, ph) in enumerate(pTbins):
			for n in range(Nharmonic_diff):
				w = data['N'][nl:nh, j] * ds[exp]['flow']['N'][nl:nh] + 1e-16
				dn2 = data['Qn'][nl:nh,j,n].conjugate() * ds[exp]['flow']['Qn'][nl:nh, n] / w
				avg_dn2, std_avg_dn2 = weighted_mean_std(dn2, w)
				vn[i, j, n] = avg_dn2/vnref['obs'][i,n]
				vn_err[i, j, n] = std_avg_dn2/vnref['obs'][i,n]
	return {'Name': 'vn2', 'cenM': cenM, 'pTM' : pTM,
			'obs': vn, 'err': vn_err}

def load_and_compute(inputfile, plot=False):
	res = np.fromfile(inputfile, dtype=result_dtype)
	res = np.array(sorted(res, key=lambda x: x['ALICE']['dNch_deta'], reverse=True))
	
	# dNdeta
	cenb = np.array([[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]])
	info = calculate_dNdeta(res, 'ALICE', cenb)
	if plot:
		plt.subplot(2,4,1)
		plt.errorbar(info['cenM'], info['obs'], yerr=info['err'], fmt='ro-')
		plt.ylim(0,2000)
		plt.xlabel(r'Centrality (%)', fontsize=7)
		plt.ylabel(r'charged $dN/d\eta$', fontsize=7)
		
	# dETdeta
	cenb = np.array([[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]])
	info = calculate_dETdeta(res, 'ALICE', cenb)
	if plot:
		plt.subplot(2,4,2)
		plt.errorbar(info['cenM'], info['obs'], yerr=info['err'], fmt='bo-')
		plt.ylim(0,2000)
		plt.xlabel(r'Centrality (%)', fontsize=7)
		plt.ylabel(r'$dE_T/d\eta$', fontsize=7)
		
	# dN(pid)/deta
	cenb = np.array([[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]])
	info = calculate_dNdy(res, 'ALICE', cenb)
	if plot:
		plt.subplot(2,4,3)
		for (s, _) in pi_K_p:
			plt.errorbar(info['cenM'], info['obs'][s], yerr=info['err'][s], fmt='o-', label=s)
		plt.ylim(0,2000)
		plt.xlabel(r'Centrality (%)', fontsize=7)
		plt.ylabel(r'$dN/dy$', fontsize=7)
		plt.legend()
		
	# mean-pT
	cenb = np.array([[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]])
	info = calculate_mean_pT(res, 'ALICE', cenb)
	if plot:
		plt.subplot(2,4,4)
		for (s, _) in species:
			plt.errorbar(info['cenM'], info['obs'][s], yerr=info['err'][s], fmt='o-', label=s)
		plt.ylim(0,2)
		plt.xlabel(r'Centrality (%)', fontsize=7)
		plt.ylabel(r'$\langle p_T\rangle$', fontsize=7)
		plt.legend()
		
	# mean-pT-fluct
	cenb = np.array([[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]])
	info = calculate_mean_pT_fluct(res, 'ALICE', cenb)
	if plot:
		plt.subplot(2,4,5)
		plt.errorbar(info['cenM'], info['obs'], yerr=info['err'], fmt='o-')
		plt.ylim(0,2)
		plt.xlabel(r'Centrality (%)', fontsize=7)
		plt.ylabel(r'$\delta p_T/p_T$', fontsize=7)
		plt.legend()

	# vn
	cenb = np.array([[0,5],[5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80],[80,90]])
	info = calculate_vn(res, 'ALICE', cenb)
	if plot:
		plt.subplot(2,4,6)
		for n in range(Nharmonic):
			plt.errorbar(info['cenM'], info['obs'][:,n], yerr=info['err'][:,n], fmt='o-', label=r"$v_{:d}$".format(n+1))
		plt.ylim(0,.15)
		plt.xlabel(r'Centrality (%)', fontsize=7)
		plt.ylabel(r'charged $dN/d\eta$', fontsize=7)
		plt.legend()
		
	# vn-diff
	cenb = np.array([[30,40]])
	pTbins = np.array([[0,.2], [.2,.4], [.4,.6],[.6,.8],[.8,1.],
			[1.,1.2], [1.2,1.5], [1.5,2.], [2.,2.5], [2.5,3]])
	info = calculate_diff_vn(res, 'ALICE', cenb, pTbins, pid='chg')
	if plot:
		plt.subplot(2,4,7)
		for n in range(Nharmonic_diff):
			plt.errorbar(info['pTM'], info['obs'][0,:,n], yerr=info['err'][0,:,n], fmt='o-', label=r"$v_{:d}$".format(n+1))
		plt.ylim(0,.15)
		plt.xlabel(r'$p_T$ [GeV]', fontsize=7)
		plt.ylabel(r'charged $v_{:d}$'.format(n+1), fontsize=7)
		plt.legend()
	
	if plot:
		plt.tight_layout(True)
		plt.show()

if __name__ == '__main__':	
	for file in glob.glob(sys.argv[1]):
		fig = plt.figure(figsize=(8,6))
		load_and_compute(file, plot=True)



