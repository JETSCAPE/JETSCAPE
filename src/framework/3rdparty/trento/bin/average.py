#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Sat 09 Jul 2016 12:30:10 AM CEST

import matplotlib.pyplot as plt
import numpy as np
import h5py

def weighted_average(h5, cent, n=2, kind='ecc'):
    '''Calc the weighted average for bulk or sd
    Params:
        :param h5: h5py file
        :param cent: centrality bin, example cent="20_30"
        :param n: harmonic order, n=2 or 3; default=2
        :param kind: string, kind='ecc', 'sd' or 'dNdEta' or 'dNdPt'

    Return:
        Average using weight from the probability of this eccentricity bin
    '''
    path = '%s/ecc%s/'%(cent, n)

    # after the 19th has tiny or no probability
    # only use the first 18 eccentricity bins;
    probability = h5[path + 'probability'][...][:18]
    probability = probability / probability.sum()

    if kind == 'ecc':
        # calc the mean eccentricity for the first 18 bins
        eccn_edges = h5[path + 'bin_edges'][...][:19]
        eccn = 0.5 * (eccn_edges[1:] + eccn_edges[:-1])
        return (probability * eccn).sum()
    elif kind == 'sd':
        nx, ny = 100, 100
        sd_18bins = np.empty((18, nx, ny))
        for i in range(18):
            sd_18bins[i] = h5[path + 'sd_bin%s'%i][...]

        return np.average(sd_18bins, axis=0, weights=probability)
    elif kind == 'dndeta':
        neta = 50
        dndeta = np.empty((18, neta, 2))
        for i in range(18):
            dndeta[i] = h5[path + 'spec_bin%s/dndeta'%i][...]

        eta_dndeta = np.average(dndeta, axis=0, weights=probability)
        return eta_dndeta[:, 0], eta_dndeta[:, 1]
        




if __name__=='__main__':
    with h5py.File('pbpb2p76.h5', 'r') as h5:
        #ecc2 = weighted_average(h5, cent='30_40', n=2, kind='ecc')
        #print(ecc2)

        #sd = weighted_average(h5, cent='30_40', n=2, kind='sd')
        #plt.imshow(sd)
        #plt.show()

        eta, dndeta_30_40 = weighted_average(h5, cent='30_40', n=2, kind='dndeta')
        eta, dndeta_0_10 = weighted_average(h5, cent='0_10', n=2, kind='dndeta')
        eta, dndeta_10_20 = weighted_average(h5, cent='10_20', n=2, kind='dndeta')
        eta, dndeta_20_30 = weighted_average(h5, cent='20_30', n=2, kind='dndeta')

        plt.plot(eta, dndeta_0_10)
        plt.plot(eta, dndeta_10_20)
        plt.plot(eta, dndeta_20_30)
        plt.plot(eta, dndeta_30_40)
        plt.show()

       

