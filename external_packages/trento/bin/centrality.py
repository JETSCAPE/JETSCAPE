#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Sun 05 Jun 2016 10:12:20 PM CEST

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Centrality(object):
    def __init__(self, nch_list):
        self.non_bias = np.sort(nch_list)[::-1]
        self.num_of_events = len(self.non_bias)

    def get_centrality_class(self, cent_min, cent_max):
        '''get the nch range [nch_max, nch_min] for centrality class [cent_min, cent_max]
        which can be used later to select events for this centrality bin'''
        eid = lambda cent: int(self.num_of_events * cent * 0.01)
        nch_max = self.non_bias[eid(cent_min)]
        nch_min = self.non_bias[eid(cent_max)]
        return nch_max, nch_min


def dist(nch):
    plt.hist(nch, bins=50)
    #plt.yscale('log', nonposy='clip')

def e2_e3_scatter(ecc2, ecc3):
    plt.scatter(ecc2, ecc3)
    plt.show()

def Nch_vs_Npart(Nch, Npart):
    plt.scatter(Npart, Nch/(0.5*Npart))
    plt.show()



if __name__=='__main__':
    dat = np.loadtxt('pbpb_1million.txt')
    #dat = pd.read_csv('pbpb_1million.txt', dtype=np.float32, sep=' ', header=None).values

    entropy = dat[:, 3]

    #dist(entropy)

    centrality = Centrality(entropy)
    #print('[smax, smin] for 0-5%=', centrality.get_centrality_class(0, 5))
    #print('[smax, smin] for 5-10%=', centrality.get_centrality_class(5, 10))
    #print('[smax, smin] for 90-95%=', centrality.get_centrality_class(90, 95))

    smax, smin = centrality.get_centrality_class(20, 30)

    cent_20_30_filter = (entropy < smax) & (entropy > smin)

    e2 = dat[:, 4]
    e3 = dat[:, 5]
    #e2_e3_scatter(e2, e3)

    dist(e3)




    #e2_e3_scatter(dat[:, 4], dat[:, 5])

    #npart = dat[:, 2]
    #Nch_vs_Npart(npart, entropy/(0.5*npart))


