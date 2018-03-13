#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Di 07 Jun 2016 11:39:36 CEST

import matplotlib.pyplot as plt
import numpy as np

import scipy.ndimage as matrix
from scipy.ndimage.measurements import center_of_mass

event_id = 30

log = np.loadtxt('/lustre/nyx/hyihp/lpang/trento_ini/pbpb100000.log')
event1 = log[event_id]
ixmc = event1[12]
iymc = event1[13]
phi_2 = event1[8]
phi_3 = event1[9]

ecc_2 = event1[4]
ecc_3 = event1[5]

import pandas as pd

ed = pd.read_csv('/lustre/nyx/hyihp/lpang/trento_ini/pbpb100000/%05d.dat'%event_id, sep=' ', header=None, dtype=np.float64, skiprows=13).values

print('phi_2 = ', phi_2)
print('phi_3 = ', phi_3)
print('ecc_2 = ', ecc_2)
print('ecc_3 = ', ecc_3)
print('mass center = ', ixmc, iymc)

print('center_of_mass(ed) from python=', center_of_mass(ed))


plt.subplot(131)
plt.imshow(ed, origin='lower')

plt.subplot(132)

#ed_shift = matrix.shift(ed, (-iymc+49.5, -ixmc+49.5))
#
#ed_rot = matrix.rotate(ed_shift,  phi_2 * 180. / np.pi, mode='constant', reshape=False)

from rotate import rotate

ed_rot = rotate(ed, ixmc, iymc, phi_2, nx=100, ny=100)

from ecc3_from_smotthIni import eccn

print('ecc2 after rotation is:', eccn(ed_rot, ixcm=49.5, iycm=49.5, n=2))

plt.imshow(ed_rot, origin='lower')

plt.subplot(133)

#ed_rot3 = matrix.rotate(ed_shift, phi_3 * 180. / np.pi, mode='constant', reshape=False)

ed_rot3 = rotate(ed, ixmc, iymc, phi_3)

print('ecc3 after rotation is:', eccn(ed_rot3, ixcm=49.5, iycm=49.5, n=3))
plt.imshow(ed_rot3, origin='lower')

plt.show()
