import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import math
import sys

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
whichSampler = sys.argv[3]

#the range of psedoraprapidity cut
etamax = 1.0

#the range of rapidity cut
ymax = 0.5

for sample in range(0, nsamples):
    hadron_file = sample_dir + 'sample_' + str(sample) + '/final' + whichSampler + 'Hadrons.dat'
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


#pT bins
delta_pT = 0.2
pT_bins_eqwidth = []
for j in range(0, 20):
    pT_bins_eqwidth.append( j * delta_pT)
np.savetxt("plots/pT_bins.dat", pT_bins_eqwidth)

#phi bins
phi_bins = []
n_phi_bins = 20
delta_phi = 2.0 * math.pi / n_phi_bins
for j in range(0, n_phi_bins):
    phi_bins.append( j * delta_phi )
np.savetxt("plots/phi_bins.dat", phi_bins)

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
np.savetxt("plots/t_bins.dat", t_bins)

#x bins
delta_x = 1.0
x_bins = []
for j in range(-10, 10):
    x_bins.append( j * delta_x)
np.savetxt("plots/x_bins.dat", x_bins)

#y bins
delta_y = 0.4
y_bins = []
for j in range(-15, 15):
    y_bins.append( j * delta_y)
np.savetxt("plots/y_bins.dat", y_bins)

#plot pion
nPi, binsPi, patchesPi = plt.hist(pi_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Pion dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots/211_spectra.pdf')
#plt.show()
plt.close()
np.savetxt("plots/pi_dN_dpT.dat", nPi)

#plot kaon
nK, binsK, patchesK = plt.hist(k_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Kaon dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots/321_spectra.pdf')
#plt.show()
plt.close()
np.savetxt("plots/k_dN_dpT.dat", nK)

#plot proton
nP, binsP, patchesP = plt.hist(p_pT_mid, bins = pT_bins_eqwidth, normed=False)
plt.title("Proton dN/dpT")
plt.xlabel("pT (GeV)")
plt.savefig('plots/2212_spectra.pdf')
#plt.show()
plt.close()
np.savetxt("plots/p_dN_dpT.dat", nP)

#examine dN/dphi
n_pi_phi, bins_pi_phi, patches_pi_phi = plt.hist(pi_phi_mid, bins = phi_bins, normed=False)
plt.title("Pion dN/dphi")
plt.xlabel("phi")
plt.savefig('plots/211_dNdphi.pdf')
#plt.show()
plt.close()
np.savetxt("plots/pi_dN_dphi.dat", n_pi_phi)

n_k_phi, bins_k_phi, patches_k_phi = plt.hist(k_phi_mid, bins = phi_bins, normed=False)
plt.title("Kaon dN/dphi")
plt.xlabel("phi")
plt.savefig('plots/321_dNdphi.pdf')
#plt.show()
plt.close()
np.savetxt("plots/k_dN_dphi.dat", n_k_phi)

n_p_phi, bins_p_phi, patches_p_phi = plt.hist(p_phi_mid, bins = phi_bins, normed=False)
plt.title("Proton dN/dphi")
plt.xlabel("phi")
plt.savefig('plots/2212_dNdphi.pdf')
#plt.show()
plt.close()
np.savetxt("plots/p_dN_dphi.dat", n_p_phi)

#examine dN/deta
plt.hist(pi_eta, bins = eta_bins, normed=False)
plt.title("Pion dN/deta")
plt.xlabel("eta")
plt.savefig('plots/211_dNdeta.pdf')
#plt.show()
plt.close()

plt.hist(k_eta, bins = eta_bins, normed=False)
plt.title("Kaon dN/deta")
plt.xlabel("eta")
plt.savefig('plots/321_dNdeta.pdf')
#plt.show()
plt.close()

plt.hist(p_eta, bins = eta_bins, normed=False)
plt.title("Proton dN/deta")
plt.xlabel("eta")
plt.savefig('plots/2212_dNdeta.pdf')
#plt.show()
plt.close()

#examine dN/dy
n_pi_y, bins_pi_y, patches_pi_y = plt.hist(pi_y, bins = y_bins, normed=False)
plt.title("Pion dN/dy")
plt.xlabel("y")
plt.savefig('plots/211_dNdy.pdf')
#plt.show()
plt.close()
np.savetxt("plots/211_dN_dy.dat", n_pi_y)

plt.hist(k_y, bins = y_bins, normed=False)
plt.title("Kaon dN/dy")
plt.xlabel("y")
plt.savefig('plots/321_dNdy.pdf')
#plt.show()
plt.close()

plt.hist(p_y, bins = y_bins, normed=False)
plt.title("Proton dN/dy")
plt.xlabel("y")
plt.savefig('plots/2212_dNdy.pdf')
#plt.show()
plt.close()

#examine dN/dt
n_t, bins_t, patches_t = plt.hist(all_t_mid, bins = t_bins, normed=False)
plt.title("dN / dt")
plt.xlabel("t")
plt.savefig('plots/dNdt.pdf')
#plt.show()
plt.close()
np.savetxt("plots/dN_dt.dat", n_t)

#examine dN/dx
n_x, bins_x, patches_x = plt.hist(all_x_mid, bins = x_bins, normed=False)
plt.title("dN / dx")
plt.xlabel("x")
plt.savefig('plots/dNdx.pdf')
#plt.show()
plt.close()
np.savetxt("plots/dN_dx.dat", n_x)
