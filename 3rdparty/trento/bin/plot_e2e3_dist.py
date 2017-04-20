#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Mon 06 Jun 2016 12:17:35 CEST

from matplotlib.pyplot import *
import numpy as np


from e2_e3_dist import *

dat = np.loadtxt('pbpb_1million.txt')
entropy = dat[:, 3]
centrality = Centrality(entropy)
smax, smin = centrality.get_centrality_class(20, 30)
cent_20_30_filter = (entropy < smax) & (entropy > smin)
c20_30 = dat[cent_20_30_filter]
e2, e3 = c20_30[:, 4], c20_30[:, 5]

subplot(121)
dist(e2)
xlabel(r'$\epsilon_2$')
ylabel(r'$N_{event}$')

subplot(122)
dist(e3)
xlabel(r'$\epsilon_3$')

from common_plotting import smash_style

smash_style.set()

suptitle(r'$Pb+Pb\ \sqrt{s_{NN}}=2.76\ TeV\ 20-30\%$', size=40)

plt.show()
