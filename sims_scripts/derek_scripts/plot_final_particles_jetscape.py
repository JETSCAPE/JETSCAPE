import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import sys

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

#midrapidity
pi_pT_mid = []
k_pT_mid = []
p_pT_mid = []


particle_list = pd.read_csv(sys.argv[1], sep=' ')
print("Reading file : " + sys.argv[1] )
mcid = particle_list['mcid']
E  = particle_list['energy']
pT = particle_list['pT']
phi = particle_list['phi']
eta = particle_list['eta']

#3d momentum space lists
for i in range(0, len(E) ):
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
etamax = 10.0
for i in range(0, len(pi_pT) ):
    if ( abs(pi_eta[i]) < etamax ):
        pi_pT_mid.append( pi_pT[i] )
for i in range(0, len(k_pT) ):
    if ( abs(k_eta[i]) < etamax ):
        k_pT_mid.append( k_pT[i] )
for i in range(0, len(p_pT) ):
    if ( abs(p_eta[i]) < etamax ):
        p_pT_mid.append( p_pT[i] )

#pT bins
delta_pT = 0.2
pT_bins_eqwidth = []
for j in range(0, 20):
    pT_bins_eqwidth.append( j * delta_pT)

#plot pion
nPi, binsPi, patchesPi = plt.hist(pi_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Pion dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots/211_spectra.pdf')
plt.show()

#plot kaon
nK, binsK, patchesK = plt.hist(k_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Kaon dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots/321_spectra.pdf')
plt.show()

#plot proton
nP, binsP, patchesP = plt.hist(p_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Proton dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots/2212_spectra.pdf')
plt.show()

#examine dN/dphi
plt.hist(pi_phi, bins = 'auto', normed=False)
plt.title("Pion 1/N_part dN/dphi")
plt.xlabel("phi")
plt.savefig('plots/211_dNdphi.pdf')
plt.show()

plt.hist(k_phi, bins = 'auto', normed=False)
plt.title("Kaon 1/N_part dN/dphi")
plt.xlabel("phi")
plt.savefig('plots/321_dNdphi.pdf')
plt.show()

plt.hist(p_phi, bins = 'auto', normed=False)
plt.title("Proton 1/N_part dN/dphi")
plt.xlabel("phi")
plt.savefig('plots/2212_dNdphi.pdf')
plt.show()

#examine dN/deta
plt.hist(pi_eta, bins = 'auto', normed=False)
plt.title("Pion dN/deta")
plt.xlabel("eta")
plt.savefig('plots/211_dNdeta.pdf')
plt.show()

plt.hist(k_eta, bins = 'auto', normed=False)
plt.title("Kaon dN/deta")
plt.xlabel("eta")
plt.savefig('plots/321_dNdeta.pdf')
plt.show()

plt.hist(p_eta, bins = 'auto', normed=False)
plt.title("Proton dN/deta")
plt.xlabel("eta")
plt.savefig('plots/2212_dNdeta.pdf')
plt.show()
