#/usr/bin/env python
#author: lgpang
#email: lgpang@qq.com
#createTime: Mon 06 Jun 2016 13:58:16 CEST

import matplotlib.pyplot as plt
import numpy as np
import h5py
from centrality import Centrality
from rotate import rotate


#path = '/lustre/nyx/hyihp/lpang/trento_ini/auau1million/'

#path = '/lustre/nyx/hyihp/lpang/trento_ini/pbpb5020/'
#dat = np.loadtxt('/lustre/nyx/hyihp/lpang/trento_ini/pbpb5020.log')

#collision_system = 'pbpb1million'
collision_system = 'pbpb5020'

path = '/lustre/nyx/hyihp/lpang/trento_ini/%s'%collision_system

dat = np.loadtxt('/lustre/nyx/hyihp/lpang/trento_ini/%s.log'%collision_system)

entropy = dat[:, 3]

# construct the centrality from nch/entropy distribution
centrality = Centrality(entropy)


def add_e2e3(h5, cent_min, cent_max):
    '''add ecc2, ecc3 distribution for this centrlaity bin in hdf5 file
       add smin, smax for this centrality bin in hdf5 file'''
    smax, smin = centrality.get_centrality_class(cent_min, cent_max)
    filter = (entropy < smax) & (entropy > smin)
    c20_30 = dat[filter]
    e2, e3 = c20_30[:, 4], c20_30[:, 5]

    e2_hist, bins_edge2 = np.histogram(e2, bins=20, density=True)
    h5.create_dataset('%s_%s/ecc2/bin_edges'%(cent_min, cent_max), data=bins_edge2)
    h5.create_dataset('%s_%s/ecc2/probability'%(cent_min, cent_max), data=e2_hist)
    h5['%s_%s'%(cent_min, cent_max)].attrs['smin'] = smin
    h5['%s_%s'%(cent_min, cent_max)].attrs['smax'] = smax
    
    e3_hist, bins_edge3 = np.histogram(e3, bins=20, density=True)
    h5.create_dataset('%s_%s/ecc3/bin_edges'%(cent_min, cent_max), data=bins_edge3)
    h5.create_dataset('%s_%s/ecc3/probability'%(cent_min, cent_max), data=e3_hist)


def add_avg_ed(h5, cent="20_30", n=2):
    '''add weighted average sd to h5, in the eccn/sd_mean'''
    from average import weighted_average
    sd_mean = weighted_average(h5, cent, n, kind='sd')
    ecc_mean = weighted_average(h5, cent, n, kind='ecc')
    ed_path = cent + '/ecc%s/sd_mean'%n
    h5.create_dataset(ed_path, data = sd_mean)
    h5[ed_path].attrs['ecc_mean'] = ecc_mean
    print("<ecc%s>="%n)



def add_ed(h5, cent_min, cent_max, nx=100, ny=100):
    ''' add the <ed>_2 which is the ed distribution rotated by participant plane phi_2
    and add the <ed>_3 which is the ed distribution rotated by participant plane phi_3
    for events in this centrality class, for all the ecc2, ecc3 bins in their probability
    distribution function'''

    add_e2e3(h5, cent_min, cent_max)
    
    cent = '%s_%s'%(cent_min, cent_max)
    smin = h5[cent].attrs['smin']
    smax = h5[cent].attrs['smax']

    filter1 = (entropy < smax) & (entropy > smin)
    events = dat[filter1]

    psi_2 = events[:, 8]
    psi_3 = events[:, 9]
    ixmc = events[:, 12]
    iymc = events[:, 13]

    event_id = events[:, 0]

    ecc2, ecc3 = events[:, 4], events[:, 5]

    e2_edges = h5[cent + '/ecc2/bin_edges'][...]
    e3_edges = h5[cent + '/ecc3/bin_edges'][...]
    e2_probs = h5[cent + '/ecc2/probability'][...]
    e3_probs = h5[cent + '/ecc3/probability'][...]

    from ecc3_from_smotthIni import eccn

    # add avg ini condition for one ecc2 bin
    #for idx, prob in enumerate(e2_probs):
    for idx in range(len(e2_probs)-1):
        ed_tot = np.zeros(shape=(ny, nx))
        ecc_min = e2_edges[idx]
        ecc_max = e2_edges[idx + 1]
        filter2 = (ecc2 >= ecc_min) & (ecc2 < ecc_max)
        id2 = event_id[filter2]
        for idx2, eid in enumerate(id2):
            ed_old = np.loadtxt('%s/%06d.dat'%(path, eid))
            ed_new = rotate(ed_old, ixmc[filter2][idx2], iymc[filter2][idx2], psi_2[filter2][idx2])
            if ecc2[filter2][idx2] < ecc_min or ecc2[filter2][idx2] > ecc_max:
                print("filter2 error")

            ed_tot += ed_new

        num_of_events_in_eccbin = float(len(id2))
        ed_path = cent + '/ecc2/sd_bin%d'%idx
        h5.create_dataset(ed_path, data = ed_tot / num_of_events_in_eccbin)
        h5[ed_path].attrs['num_of_events'] = num_of_events_in_eccbin
        print('ecc2 bin ', idx, ' finished, num of events in this bin=', num_of_events_in_eccbin)

    # add avg ini condition for one ecc3 bin
    for idx, prob in enumerate(e3_probs):
        ed_tot = np.zeros(shape=(nx, ny))
        ecc_min = e3_edges[idx]
        ecc_max = e3_edges[idx + 1]
        filter3 = (ecc3 >= ecc_min) & (ecc3 < ecc_max)
        id2 = event_id[filter3]
        for idx2, eid in enumerate(id2):
            ed_old = np.loadtxt('%s/%06d.dat'%(path, eid))
            ed_tot += rotate(ed_old, ixmc[filter3][idx2], iymc[filter3][idx2], psi_3[filter3][idx2])

        num_of_events_in_eccbin = float(len(id2))
        ed_path = cent + '/ecc3/sd_bin%d'%idx
        h5.create_dataset(ed_path, data = ed_tot / num_of_events_in_eccbin)
        h5[ed_path].attrs['num_of_events'] = num_of_events_in_eccbin
        print('ecc3 bin ', idx, ' finished, num of events in this bin=', num_of_events_in_eccbin)

    # do weighted average for <sd> and <ecc> according to their probability
    add_avg_ed(h5, cent, n=2)
    add_avg_ed(h5, cent, n=3)


with h5py.File('%s_new.h5'%collision_system, 'w') as h5:
    h5.create_dataset('events_info', data=dat)

    cents = [(0, 5), (0, 10), (5, 10), (10, 20), (20, 30),
             (0, 20), (10, 30), (30, 40), (30, 50), (40, 50), (0, 80)]

    for cent_min, cent_max in cents:
        add_ed(h5, cent_min, cent_max)


