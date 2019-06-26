#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import h5py, sys

def help():
	
	"""
  This script plots 
    (1) transverse intergated entropy/density as function of rapidity
    (2) Projection of entropy/density onto x-y (eta=0) plane
    (3) ............................. onto y-eta (x=0) plane.
    (4) ............................. onto x-eta (y=0) plane.

  Usage:  
    {:s} trento3d-hdf5-output [list-of-event-id-to-convert]

  For example, to view all events:
    {:s} ic.hdf5
  To view only events #2 and #3:
    {:s} ic.hdf5 2 3
"""
	print(help.__doc__.format(__file__, __file__, __file__))

def plot(dataset):
	fig, axes = plt.subplots(
		nrows=2, ncols=2,
		figsize=(8, 7)
	)
	ax = axes[0,0]
	field = dataset.value
	Nx, Ny, Neta = field.shape
	deta = dataset.attrs['deta']
	dxy = dataset.attrs['dxy']
	Lx, Ly, Leta = Nx*dxy/2., Ny*dxy/2., (Neta-1)*deta/2.
	x = np.linspace(-Lx, Lx, Nx)
	y = np.linspace(-Ly, Ly, Ny)
	eta = np.linspace(-Leta, Leta, Neta)
	ax.plot(eta, np.sum(field, axis=(0,1))*dxy*dxy)
	ax.set_xlabel(r'$\eta$')
	ax.set_ylabel(r'$dS/d\eta$')
	ax.set_title('Tranverse integrated entropy')
	Nxmid, Nymid, Netamid = int((Nx-1)/2), int((Ny-1)/2), int((Neta-1)/2)
	for ax, projection, ranges, xlabel, ylabel, title in zip(
		axes.flatten()[1:], 
		[field[:,:,Netamid], field[Nxmid,:,:], field[:,Nymid,:]],
		[[-Lx, Lx, -Ly, Ly], [-Ly, Ly, -Leta, Leta], [-Lx, Lx, -Leta, Leta]],
		[r'$y$ [fm]', r'$\eta$', r'$\eta$'], 
		[r'$x$ [fm]', r'$y$ [fm]', r'$x$ [fm]'],
		[r'$x$-$y$ projection, $\eta=0$',r'$y$-$\eta$ projection, $x=0$',
		 r'$x$-$\eta$ projection, $y=0$']
		):
		ax.contourf(projection, extent=ranges)
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_title(title)
		ax.axis('equal')
	plt.subplots_adjust(wspace=0.4, hspace=0.4)
	plt.show()

def main():
	if len(sys.argv) <= 1:
		help()
		exit()
	f = h5py.File(sys.argv[1], 'r')
	elist = ['event_{}'.format(index) for index in sys.argv[2:]] \
			if len(sys.argv)>=3 else list(f.keys())
	for eid in elist:
		plot(f[eid])

if __name__ == '__main__':
	main()
