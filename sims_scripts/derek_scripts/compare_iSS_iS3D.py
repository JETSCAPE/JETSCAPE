import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import math
import sys

#iSS
#######################
#species dependent info
#pion 211
pi_pT  = []
pi_eta = []
pi_phi = []

#kaon 321
k_pT  = []
k_eta = []
k_phi = []

#proton 2212
p_pT  = []
p_eta = []
p_phi = []

#after cut on pseudorap
pi_pT_mid = []
k_pT_mid = []
p_pT_mid = []

pi_phi_mid = []
k_phi_mid = []
p_phi_mid = []

t_mid = []
x_mid = []

iSS_sample_dir = sys.argv[1]
nsamples = int(sys.argv[3])

for sample in range(0, nsamples):
    hadron_file = iSS_sample_dir + 'sample_' + str(sample) + '/finaliSSHadrons.dat'
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
        elif ( mcid[i] == 321 ):
            k_pT.append( pT[i] )
            k_eta.append( eta[i] )
            k_phi.append( phi[i] )
        elif ( mcid[i] == 2212 ):
            p_pT.append( pT[i] )
            p_eta.append( eta[i] )
            p_phi.append( phi[i] )

#the range of psedoraprapidity cut
etamax = 3.0
for i in range(0, len(pi_pT) ):
    if ( abs(pi_eta[i]) < etamax ):
        pi_pT_mid.append( pi_pT[i] )
        pi_phi_mid.append( pi_phi[i] )
for i in range(0, len(k_pT) ):
    if ( abs(k_eta[i]) < etamax ):
        k_pT_mid.append( k_pT[i] )
        k_phi_mid.append( k_phi[i] )
for i in range(0, len(p_pT) ):
    if ( abs(p_eta[i]) < etamax ):
        p_pT_mid.append( p_pT[i] )
        p_phi_mid.append( p_phi[i] )

for i in range(0, len(eta) ):
    if ( abs(eta[i]) < etamax ):
        t_mid.append( t[i] )
        x_mid.append( x[i] )

#pT bins
delta_pT = 0.2
pT_bins_eqwidth = []
for j in range(0, 20):
    pT_bins_eqwidth.append( j * delta_pT)

#phi bins
phi_bins = []
n_phi_bins = 20
delta_phi = 2.0 * math.pi / n_phi_bins
for j in range(0, n_phi_bins):
    phi_bins.append( j * delta_phi )

#eta bins
delta_eta = 0.2
eta_bins = []
for j in range(-10, 10):
    eta_bins.append( j * delta_eta)

#time bins
delta_t = 0.2
t_bins = []
for j in range(-10, 10):
    t_bins.append( j * delta_t)

#x bins
delta_x = 1.0
x_bins = []
for j in range(-10, 10):
    x_bins.append( j * delta_x)

#plot pion
nPi, binsPi, patchesPi = plt.hist(pi_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Pion dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots_iSS/211_spectra.pdf')
#plt.show()
plt.close()

#plot kaon
nK, binsK, patchesK = plt.hist(k_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Kaon dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots_iSS/321_spectra.pdf')
#plt.show()
plt.close()

#plot proton
nP, binsP, patchesP = plt.hist(p_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Proton dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots_iSS/2212_spectra.pdf')
#plt.show()
plt.close()

#examine dN/dphi
plt.hist(pi_phi_mid, bins = phi_bins, normed=False)
plt.title("Pion 1/N_part dN/dphi")
plt.xlabel("phi")
plt.savefig('plots_iSS/211_dNdphi.pdf')
#plt.show()
plt.close()

plt.hist(k_phi_mid, bins = phi_bins, normed=False)
plt.title("Kaon 1/N_part dN/dphi")
plt.xlabel("phi")
plt.savefig('plots_iSS/321_dNdphi.pdf')
#plt.show()
plt.close()

plt.hist(p_phi_mid, bins = phi_bins, normed=False)
plt.title("Proton 1/N_part dN/dphi")
plt.xlabel("phi")
plt.savefig('plots_iSS/2212_dNdphi.pdf')
#plt.show()
plt.close()

#examine dN/deta
plt.hist(pi_eta, bins = eta_bins, normed=False)
plt.title("Pion dN/deta")
plt.xlabel("eta")
plt.savefig('plots_iSS/211_dNdeta.pdf')
#plt.show()
plt.close()

plt.hist(k_eta, bins = eta_bins, normed=False)
plt.title("Kaon dN/deta")
plt.xlabel("eta")
plt.savefig('plots_iSS/321_dNdeta.pdf')
#plt.show()
plt.close()

plt.hist(p_eta, bins = eta_bins, normed=False)
plt.title("Proton dN/deta")
plt.xlabel("eta")
plt.savefig('plots_iSS/2212_dNdeta.pdf')
#plt.show()
plt.close()

#examine dN/dt
plt.hist(t_mid, bins = t_bins, normed=False)
plt.title("dN / dt")
plt.xlabel("t")
plt.savefig('plots_iSS/dNdt.pdf')
#plt.show()
plt.close()

#examine dN/dx
plt.hist(x_mid, bins = x_bins, normed=False)
plt.title("dN / dx")
plt.xlabel("x")
plt.savefig('plots_iSS/dNdx.pdf')
#plt.show()
plt.close()
