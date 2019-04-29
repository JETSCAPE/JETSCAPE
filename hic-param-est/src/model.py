"""
Computes model observables to match experimental data.
Prints all model data when run as a script.

Model data files are expected with the file structure
:file:`model_output/{design}/{system}/{design_point}.dat`, where
:file:`{design}` is a design type, :file:`{system}` is a system string, and
:file:`{design_point}` is a design point name.

For example, the structure of my :file:`model_output` directory is ::

    model_output
    ├── main
    │   ├── PbPb2760
    │   │   ├── 000.dat
    │   │   └── 001.dat
    │   └── PbPb5020
    │       ├── 000.dat
    │       └── 001.dat
    └── validation
        ├── PbPb2760
        │   ├── 000.dat
        │   └── 001.dat
        └── PbPb5020
            ├── 000.dat
            └── 001.dat

I have two design types (main and validation), two systems, and my design
points are numbered 000-499 (most numbers omitted for brevity).

Data files are expected to have the binary format created by my `heavy-ion
collision event generator
<https://github.com/jbernhard/heavy-ion-collisions-osg>`_.

Of course, if you have a different data organization scheme and/or format,
that's fine.  Modify the code for your needs.
"""

import logging
from pathlib import Path
import pickle

from hic import flow
import numpy as np
from sklearn.externals import joblib

from . import workdir, cachedir, systems, lazydict, expt
from .design import Design


def pT_fluct(events):
    """
    Compute the relative mean pT fluctuation \sqrt{C_m}/M(p_T)_m, defined in
    Eqs. (2-3) of the ALICE paper https://inspirehep.net/record/1307102.

    """
    N, sum_pT, sum_pTsq = [
        events['pT_fluct'][k] for k in
        ['N', 'sum_pT', 'sum_pTsq']
    ]

    Npairs = .5*N*(N - 1)
    M = sum_pT.sum() / N.sum()

    # This is equivalent to the sum over pairs in Eq. (2).  It may be derived
    # by using that, in general,
    #
    #   \sum_{i,j>i} a_i a_j = 1/2 [(\sum_{i} a_i)^2 - \sum_{i} a_i^2].
    #
    # That is, the sum over pairs (a_i, a_j) may be re-expressed in terms of
    # the sum of a_i and sum of squares a_i^2.  Applying this to Eq. (2) and
    # collecting terms yields the following expression.
    C = (
        .5*(sum_pT**2 - sum_pTsq) - M*(N - 1)*sum_pT + M**2*Npairs
    ).sum() / Npairs.sum()

    return np.sqrt(C)/M


# TODO move this symmetric cumulant code to hic

def csq(x):
    """
    Return the absolute square |x|^2 of a complex array.

    """
    return (x*x.conj()).real


def corr2(Qn, N):
    """
    Compute the two-particle correlation <v_n^2>.

    """
    return (csq(Qn) - N).sum() / (N*(N - 1)).sum()


def symmetric_cumulant(events, m, n, normalize=False):
    """
    Compute the symmetric cumulant SC(m, n).

    """
    N = np.asarray(events['flow']['N'], dtype=float)
    Q = dict(enumerate(events['flow']['Qn'].T, start=1))

    cm2n2 = (
        csq(Q[m]) * csq(Q[n])
        - 2*(Q[m+n] * Q[m].conj() * Q[n].conj()).real
        - 2*(Q[m] * Q[m-n].conj() * Q[n].conj()).real
        + csq(Q[m+n]) + csq(Q[m-n])
        - (N - 4)*(csq(Q[m]) + csq(Q[n]))
        + N*(N - 6)
    ).sum() / (N*(N - 1)*(N - 2)*(N - 3)).sum()

    cm2 = corr2(Q[m], N)
    cn2 = corr2(Q[n], N)

    sc = cm2n2 - cm2*cn2

    if normalize:
        sc /= cm2*cn2

    return sc


# fully specify numeric data types, including endianness and size, to
# ensure consistency across all machines
float_t = '<f8'
int_t = '<i8'
complex_t = '<c16'


class ModelData:
    """
    Helper class for event-by-event model data.  Reads binary data files and
    computes centrality-binned observables.

    """
    species = ['pion', 'kaon', 'proton', 'Lambda', 'Sigma0', 'Xi', 'Omega']

    #: The expected binary data type.
    dtype = np.dtype([
        ('initial_entropy', float_t),
        ('nsamples', int_t),
        ('dNch_deta', float_t),
        ('dET_deta', float_t),
        ('dN_dy', [(s, float_t) for s in species]),
        ('mean_pT', [(s, float_t) for s in species]),
        ('pT_fluct', [('N', int_t), ('sum_pT', float_t), ('sum_pTsq', float_t)]),
        ('flow', [('N', int_t), ('Qn', complex_t, 8)]),
    ])

    def __init__(self, *files):
        # read each file using the above dtype and treat each as a minimum-bias
        # event sample
        def load_events(f):
            logging.debug('loading %s', f)
            d = np.fromfile(str(f), dtype=self.dtype)
            d.sort(order='dNch_deta')
            return d

        self.events = [load_events(f) for f in files]

    def observables_like(self, data, *keys):
        """
        Compute the same centrality-binned observables as contained in `data`
        with the same nested dict structure.

        This function calls itself recursively, each time prepending to `keys`.

        """
        try:
            x = data['x']
            cent = data['cent']
        except KeyError:
            return {
                k: self.observables_like(v, k, *keys)
                for k, v in data.items()
            }

        def _compute_bin():
            """
            Choose a function to compute the current observable for a single
            centrality bin.

            """
            obs_stack = list(keys)
            obs = obs_stack.pop()

            if obs in ['dNch_deta', 'dET_deta']:
                return lambda events: events[obs].mean()

            if obs == 'dN_dy':
                species = obs_stack.pop()
                return lambda events: events[obs][species].mean()

            if obs == 'mean_pT':
                species = obs_stack.pop()
                return lambda events: np.average(
                    events[obs][species],
                    weights=events['dN_dy'][species]
                )

            if obs == 'pT_fluct':
                return pT_fluct

            if obs.startswith('vnk'):
                nk = obs_stack.pop()
                return lambda events: flow.Cumulant(
                    events['flow']['N'], *events['flow']['Qn'].T[1:]
                ).flow(*nk, imaginary='zero')

            if obs.startswith('sc'):
                mn = obs_stack.pop()
                return lambda events: symmetric_cumulant(
                    events, *mn, normalize='normed' in obs
                )

        compute_bin = _compute_bin()

        def compute_all_bins(events):
            n = events.size
            bins = [
                events[int((1 - b/100)*n):int((1 - a/100)*n)]
                for a, b in cent
            ]

            return list(map(compute_bin, bins))

        return dict(
            x=x, cent=cent,
            Y=np.array(list(map(compute_all_bins, self.events))).squeeze()
        )


def _data(system, dataset='main'):
    """
    Compute model observables for the given system and dataset.

    dataset may be one of:

        - 'main' (training design)
        - 'validation' (validation design)
        - 'map' (maximum a posteriori, i.e. "best-fit" point)

    """
    if dataset not in {'main', 'validation', 'map'}:
        raise ValueError('invalid dataset: {}'.format(dataset))

    files = (
        [Path(workdir, 'model_output', dataset, '{}.dat'.format(system))]
        if dataset == 'map' else
        [
            Path(workdir, 'model_output', dataset, system, '{}.dat'.format(p))
            for p in
            Design(system, validation=(dataset == 'validation')).points
        ]
    )

    cachefile = Path(cachedir, 'model', dataset, '{}.pkl'.format(system))

    if cachefile.exists():
        # use the cache unless any of the model data files are newer
        # this DOES NOT check any other logical dependencies, e.g. the
        # experimental data
        # to force recomputation, delete the cache file
        mtime = cachefile.stat().st_mtime
        if all(f.stat().st_mtime < mtime for f in files):
            logging.debug('loading observables cache file %s', cachefile)
            return joblib.load(cachefile)
        else:
            logging.debug('cache file %s is older than event data', cachefile)
    else:
        logging.debug('cache file %s does not exist', cachefile)

    logging.info(
        'loading %s/%s data and computing observables',
        system, dataset
    )

    data = expt.data[system]

    # some data are not yet available for PbPb5020
    # create dummy entries so that they are computed for the model
    if system == 'PbPb5020':
        data = dict(
            ((obs, obsdata) for obs, obsdata in
             expt.data['PbPb2760'].items() if obs not in data),
            **data
        )

    data = ModelData(*files).observables_like(data)

    logging.info('writing cache file %s', cachefile)
    cachefile.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(data, cachefile, protocol=pickle.HIGHEST_PROTOCOL)

    return data


data = lazydict(_data, 'main')
validation_data = lazydict(_data, 'validation')
map_data = lazydict(_data, 'map')


if __name__ == '__main__':
    from pprint import pprint
    for s in systems:
        d = data[s]
        print(s)
        pprint(d)
