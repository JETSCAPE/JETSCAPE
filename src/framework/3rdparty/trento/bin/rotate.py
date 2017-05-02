#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Di 07 Jun 2016 11:39:36 CEST

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.ndimage.interpolation as matrix


def rotate(ed, ixmc=50, iymc=50, phi_n=0.0, nx=100, ny=100):
    grid_cent_x = (nx - 1.) / 2.
    grid_cent_y = (ny - 1.) / 2.
    dx_dy = (grid_cent_y-iymc, grid_cent_x-ixmc)

    ed_shift = matrix.shift(ed, dx_dy)

    ed_rot = matrix.rotate(ed_shift,  - phi_n * 180. / np.pi, mode='constant', reshape=False, order=2)

    return ed_rot


if __name__ == '__main__':
    event_id = 6
    
    log = np.loadtxt('test.log')
    event1 = log[event_id]
    ixmc = event1[12]
    iymc = event1[13]
    phi_2 = event1[8]
    
    ed = pd.read_csv('test/%s.dat'%event_id, sep=' ', header=None, dtype=np.float64, skiprows=13)
    
    print('phi_2 = ', phi_2)
    print('mass center = ', ixmc, iymc)
    
    
    plt.subplot(121)
    plt.imshow(ed, origin='lower')
    
    plt.subplot(122)
    
    ed_shift = matrix.shift(ed, (-iymc+50., -ixmc+50.))
    
    ed_rot = matrix.rotate(ed_shift, - phi_2 * 180. / np.pi, mode='constant', reshape=False)
    
    plt.imshow(ed_rot, origin='lower')
    
    
    plt.show()
