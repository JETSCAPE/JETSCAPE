#!/usr/bin/env python3
import numpy as np, h5py, math
import sys, os, warnings
import matplotlib.pyplot as plt

def help():
	"""
  This script transforms the initial condition from
  (tau, x, y, eta) frame to (t, x, y, z) frame. This 
  conversion assumes longitudinal Bjorken expansion
  and transverse still.

  The entropy density in longitudinal-local frame at
  constant proper time tau0 is:
    s(tau0, x, y, eta) = s0(x,y,eta)/tau0 
                       = ns0(x,y,eta) <-- trento3d 

  Assuming Bjorken expansion for a short proper time,
    s(tau, x, y, eta) = 1/tau*s0(x,y,eta) 
                      = tau0/tau*ns0(x,y,eta) -- (1)

  Expression (1) is transformed into (t, x, y, z) via,
    s'(t0, x, y, z) = 1/tau*ns0(x, y, eta)
                tau = sqrt(t0^2 - z^2)
                eta = 1/2*ln[(t0+z)/(t0-z)]
  where tau0 is absorbed into normalization of ns0
  (see --normalization option of the trento3d model).

  Usage:  
    {:s} trento3d-hdf5-output t0 [list-of-event-id-to-convert]

  For example, to convert all events to t0=1fm/c:
    {:s} ic.hdf5 1.0
  To convert only events #2 and #3 to t0=1fm/c:
    {:s} ic.hdf5 1.0 2 3
"""
	print(help.__doc__.format(__file__, __file__, __file__))

tiny = 1e-9

def convert(dataset, t0):
	field = dataset.value
	Nx, Ny, Neta = field.shape
	deta = dataset.attrs['deta']
	dxy = dataset.attrs['dxy']
	Lx, Ly, Leta = Nx*dxy/2., Ny*dxy/2., (Neta-1)*deta/2.
	x = np.linspace(-Lx, Lx, Nx)
	y = np.linspace(-Ly, Ly, Ny)
	eta = np.linspace(-Leta, Leta, Neta)
	# set maximum z-range, about +/- 6-unit of rapidity
	z_max = t0*math.tanh(6.)
	zarray = np.linspace(-z_max, z_max, Neta)

	new_field = np.zeros_like(field)
	
	for iz, z in enumerate(zarray):
		eta = 0.5*np.log((t0+z)/(t0-z))
		if np.abs(eta) > Leta:
			continue
		residue, index = math.modf((eta+Leta)/deta)
		index = int(index)
		new_field[:,:,iz] = (  field[:,:,index]  *(1.-residue) \
							 + field[:,:,index+1]*residue) 	\
							 / np.sqrt(t0**2-z**2)
	return x, y, zarray, new_field

def plot(x, y, z, s):
	dx = x[1]-x[0]
	dy = y[1]-y[0]
	fig, ax = plt.subplots(1, 1)
	ax.plot(z, np.sum(s, axis=(0,1))*dx*dy, 'ro-')
	ax.set_xlabel(r'$z$ [fm]', size=15)
	ax.set_title(r'Local-rest-frame entropy in $t-z$ coordinates', size=15)
	ax.set_ylabel(r'$ds/dV$ [Arb. units.]', size=15)
	plt.subplots_adjust(wspace=0.15, bottom=0.2)
	plt.semilogy()
	plt.show()

def main():
	if len(sys.argv) <= 2:
		help()
		exit()
	f = h5py.File(sys.argv[1], 'r')
	t0 = float(sys.argv[2])
	elist = ['event_{}'.format(index) for index in sys.argv[3:]] \
			if len(sys.argv)>=4 else list(f.keys())	
	for eid in elist:
		x, y, z, s = convert(f[eid], t0)
		plot(x, y, z, s)
if __name__ == "__main__":
	main()
