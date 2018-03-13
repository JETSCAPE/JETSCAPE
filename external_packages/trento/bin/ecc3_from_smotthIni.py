#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Fr 10 Jun 2016 10:20:09 CEST

import matplotlib.pyplot as plt
import numpy as np
import h5py
import cmath
from scipy.ndimage.measurements import center_of_mass

def eccn(ed, ixcm=49.5, iycm=49.5, n=2):
    '''calc ecc_3 from shifted and rorated event average
    initial energy density distribution'''

    x = np.array(range(100)) - ixcm
    y = np.array(range(100)) - iycm

    X, Y = np.meshgrid(x, y, indexing='xy')

    print('center_of_mass=', center_of_mass(ed))

    r2 = X*X + Y*Y
    rn = np.power(r2, n/2.)

    phi = np.arctan2(Y, X)

    norm = (rn * np.exp(1j * n * phi) * ed).sum()
    denm = (rn * ed).sum()

    return abs(norm/denm), cmath.phase(norm)/float(n)
    




if __name__ == '__main__':
    n = 2
    
    with h5py.File('pbpb2p76.h5', 'r+') as h5:
        eccn_bin = h5['50_60/ecc%d/bin_edges'%n][...]
        for ibin in range(19):
            ed_avg = h5['50_60/ecc%d/ed_bin%d'%(n, ibin)][...]
            print ibin, eccn(ed_avg), eccn_bin[ibin], eccn_bin[ibin+1]
            #plt.imshow(ed_avg, origin='lower')
            #plt.show()


