#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Di 07 Jun 2016 22:44:03 CEST

import matplotlib.pyplot as plt
import numpy as np
from common_plotting import smash_style

from matplotlib.transforms import Bbox

import h5py

h5 = h5py.File('pbpb2p76.h5', 'r')

harmonic_order = 2

ecc2 = h5['20_30/ecc%s/bin_edges'%harmonic_order][...]
pro2 = h5['20_30/ecc%s/probability'%harmonic_order][...]


x = 0.5*(ecc2[:-1] + ecc2[1:])
y = pro2

fig = plt.figure(figsize=(12, 8))

ax_main = plt.gca()

plt.subplots_adjust(left=0.15, bottom=0.15)

plt.plot(x, y, 'o')
w, h = 0.2, 0.2

for i in range(len(x)-1):
    bb_data = Bbox.from_bounds(x[i]-0.5*w, y[i]-0.5*h, w, h)
    disp_coords = ax_main.transData.transform(bb_data)
    fig_coords = fig.transFigure.inverted().transform(disp_coords)
    rect = Bbox(fig_coords)

    ed = h5['20_30/ecc%s/ed_bin%d'%(harmonic_order,i)][...]

    ax = fig.add_axes(rect, frameon=False, axisbg='g', aspect=1)

    ax.patch.set_visible(False)

    ax.contourf(ed)
    ax.set_axis_off()

ax_main.set_xlabel(r'$\epsilon_%s$'%harmonic_order)
ax_main.set_ylabel(r'$f(\epsilon_%s)$'%harmonic_order)

ax_main.set_title(r'$Pb+Pb,\ \sqrt{s_{NN}}=2.76\ TeV, 20-30\%$')

smash_style.set()

plt.savefig('pbpb_20_30_ecc%s_vs_ed.pdf'%harmonic_order)

plt.show()
