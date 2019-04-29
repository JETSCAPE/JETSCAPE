"""
Downloads, processes, and stores experimental data.
Prints all data when run as a script.
"""

from collections import defaultdict
import logging
import pickle
from urllib.request import urlopen

import numpy as np
import yaml

from . import cachedir, systems


class HEPData:
    """
    Interface to a `HEPData <https://hepdata.net>`_ YAML data table.

    Downloads and caches the dataset specified by the INSPIRE record and table
    number.  The web UI for `inspire_rec` may be found at
    :file:`https://hepdata.net/record/ins{inspire_rec}`.

    If `reverse` is true, reverse the order of the data table (useful for
    tables that are given as a function of Npart).

    .. note::

        Datasets are assumed to be a function of centrality.  Other kinds of
        datasets will require code modifications.

    """
    def __init__(self, inspire_rec, table, reverse=False):
        cachefile = (
            cachedir / 'hepdata' /
            'ins{}_table{}.pkl'.format(inspire_rec, table)
        )
        name = 'record {} table {}'.format(inspire_rec, table)

        if cachefile.exists():
            logging.debug('loading from hepdata cache: %s', name)
            with cachefile.open('rb') as f:
                self._data = pickle.load(f)
        else:
            logging.debug('downloading from hepdata.net: %s', name)
            cachefile.parent.mkdir(exist_ok=True)
            with cachefile.open('wb') as f, urlopen(
                    'https://hepdata.net/download/table/'
                    'ins{}/Table{}/yaml'.format(inspire_rec, table)
            ) as u:
                self._data = yaml.load(u)
                pickle.dump(self._data, f, protocol=pickle.HIGHEST_PROTOCOL)

        if reverse:
            for v in self._data.values():
                for d in v:
                    d['values'].reverse()

    def x(self, name, case=True):
        """
        Get an independent variable ("x" data) with the given name.

        If `case` is false, perform case-insensitive matching for the name.

        """
        trans = (lambda x: x) if case else (lambda x: x.casefold())
        name = trans(name)

        for x in self._data['independent_variables']:
            if trans(x['header']['name']) == name:
                return x['values']

        raise LookupError("no x data with name '{}'".format(name))

    @property
    def cent(self):
        """
        The centrality bins as a list of (low, high) tuples.

        """
        try:
            return self._cent
        except AttributeError:
            pass

        x = self.x('centrality', case=False)

        if x is None:
            raise LookupError('no centrality data')

        try:
            cent = [(v['low'], v['high']) for v in x]
        except KeyError:
            # try to guess bins from midpoints
            mids = [v['value'] for v in x]
            width = set(a - b for a, b in zip(mids[1:], mids[:-1]))
            if len(width) > 1:
                raise RuntimeError('variable bin widths')
            d = width.pop() / 2
            cent = [(m - d, m + d) for m in mids]

        self._cent = cent

        return cent

    @cent.setter
    def cent(self, value):
        """
        Manually set centrality bins.

        """
        self._cent = value

    def y(self, name=None, **quals):
        """
        Get a dependent variable ("y" data) with the given name and qualifiers.

        """
        for y in self._data['dependent_variables']:
            if name is None or y['header']['name'].startswith(name):
                y_quals = {q['name']: q['value'] for q in y['qualifiers']}
                if all(y_quals[k] == v for k, v in quals.items()):
                    return y['values']

        raise LookupError(
            "no y data with name '{}' and qualifiers '{}'"
            .format(name, quals)
        )

    def dataset(self, name=None, maxcent=70, ignore_bins=[], **quals):
        """
        Return a dict containing:

        - **cent:** list of centrality bins
        - **x:** numpy array of centrality bin midpoints
        - **y:** numpy array of y values
        - **yerr:** subdict of numpy arrays of y errors

        `name` and `quals` are passed to `HEPData.y()`.

        Missing y values are skipped.

        Centrality bins whose upper edge is greater than `maxcent` are skipped.

        Centrality bins in `ignore_bins` [a list of (low, high) tuples] are
        skipped.

        """
        cent = []
        y = []
        yerr = defaultdict(list)

        for c, v in zip(self.cent, self.y(name, **quals)):
            # skip missing values
            # skip bins whose upper edge is greater than maxcent
            # skip explicitly ignored bins
            if v['value'] == '-' or c[1] > maxcent or c in ignore_bins:
                continue

            cent.append(c)
            y.append(v['value'])

            for err in v['errors']:
                try:
                    e = err['symerror']
                except KeyError:
                    e = err['asymerror']
                    if abs(e['plus']) != abs(e['minus']):
                        raise RuntimeError(
                            'asymmetric errors are not implemented'
                        )
                    e = abs(e['plus'])

                yerr[err.get('label', 'sum')].append(e)

        return dict(
            cent=cent,
            x=np.array([(a + b)/2 for a, b in cent]),
            y=np.array(y),
            yerr={k: np.array(v) for k, v in yerr.items()},
        )


def _data():
    """
    Curate the experimental data using the `HEPData` class and return a nested
    dict with levels

    - system
    - observable
    - subobservable
    - dataset (created by :meth:`HEPData.dataset`)

    For example, ``data['PbPb2760']['dN_dy']['pion']`` retrieves the dataset
    for pion dN/dy in Pb+Pb collisions at 2.76 TeV.

    Some observables, such as charged-particle multiplicity, don't have a
    natural subobservable, in which case the subobservable is set to `None`.

    The best way to understand the nested dict structure is to explore the
    object in an interactive Python session.

    """
    data = {s: {} for s in systems}

    # PbPb2760 and PbPb5020 dNch/deta
    for system, args, name in [
            ('PbPb2760', (880049, 1), 'D(N)/DETARAP'),
            ('PbPb5020', (1410589, 2),
             r'$\mathrm{d}N_\mathrm{ch}/\mathrm{d}\eta$'),
    ]:
        data[system]['dNch_deta'] = {None: HEPData(*args).dataset(name)}

    # PbPb2760 transverse energy
    # ignore bin 0-5 since it's redundant with 0-2.5 and 2.5-5
    dset = HEPData(1427723, 1).dataset('$E_{T}$', ignore_bins=[(0, 5)])
    dset['yerr']['sys'] = dset['yerr'].pop('sys,total')
    data['PbPb2760']['dET_deta'] = {None: dset}

    # PbPb2760 identified dN/dy and mean pT
    system = 'PbPb2760'

    for obs, table, combine_func in [
            ('dN_dy', 31, np.sum),
            ('mean_pT', 32, np.mean),
    ]:
        data[system][obs] = {}
        d = HEPData(1222333, table)
        for key, re_products in [
            ('pion', ['PI+', 'PI-']),
            ('kaon', ['K+', 'K-']),
            ('proton', ['P', 'PBAR']),
        ]:
            dsets = [
                d.dataset(RE='PB PB --> {} X'.format(i))
                for i in re_products
            ]

            data[system][obs][key] = dict(
                dsets[0],
                y=combine_func([d['y'] for d in dsets], axis=0),
                yerr={
                    e: combine_func([d['yerr'][e] for d in dsets], axis=0)
                    for e in dsets[0]['yerr']
                }
            )

    # PbPb2760 strange baryon yields
    data['PbPb2760']['dN_dy']['Lambda'] = HEPData(1243863, 23).dataset(
        RE='PB PB --> LAMBDA X'
    )

    d = HEPData(1243865, 11)
    for s in ['Xi', 'Omega']:
        data[system]['dN_dy'][s] = d.dataset(
            RE='PB PB --> ({0}- + {0}BAR+) X'.format(s.upper())
        )

    # PbPb2760 mean pT fluctuations
    d = HEPData(1307102, 6, reverse=True)
    name = r'$\sqrt{C_m}/M(p_{\rm T})_m$'
    # the table only has Npart, but they are actually 5% centrality bins
    width = 5.
    d.cent = [(n*width, (n+1)*width) for n, _ in enumerate(d.y(name))]
    data['PbPb2760']['pT_fluct'] = {None: d.dataset(name, maxcent=60)}

    # PbPb2760 and PbPb5020 flows
    for system, tables_nk in [
            ('PbPb5020', [
                (1, [(2, 2), (2, 4)]),
                (2, [(3, 2), (4, 2)]),
            ]),
            ('PbPb2760', [
                (3, [(2, 2), (2, 4)]),
                (4, [(3, 2), (4, 2)]),
            ]),
    ]:
        data[system]['vnk'] = {}

        for table, nk in tables_nk:
            d = HEPData(1419244, table)
            for n, k in nk:
                data[system]['vnk'][n, k] = d.dataset(
                    'V{}{{{}{}}}'.format(
                        n, k, ', |DELTAETA|>1' if k == 2 else ''
                    ),
                    maxcent=(70 if n == 2 else 50)
                )

    # PbPb2760 central flows vn{2}
    system, obs = 'PbPb2760', 'vnk_central'
    data[system][obs] = {}

    for n, table, sys_err_frac in [(2, 11, .025), (3, 12, .040)]:
        dset = HEPData(900651, table).dataset()
        # the (unlabeled) errors in the dataset are actually stat
        dset['yerr']['stat'] = dset['yerr'].pop('sum')
        # sys error is not provided -- use estimated fractions
        dset['yerr']['sys'] = sys_err_frac * dset['y']
        data[system][obs][n, 2] = dset

    # PbPb2760 flow correlations
    for obs, table in [
            ('sc', 1),
            ('sc_normed', 2),
            ('sc_central', 3),
            ('sc_normed_central', 4)
    ]:
        d = HEPData(1452590, table)
        data['PbPb2760'][obs] = {
            mn: d.dataset('SC({},{})'.format(*mn))
            for mn in [(3, 2), (4, 2)]
        }

    return data


#: A nested dict containing all the experimental data, created by the
#: :func:`_data` function.
data = _data()


def cov(
        system, obs1, subobs1, obs2, subobs2,
        stat_frac=1e-4, sys_corr_length=100, cross_factor=.8,
        corr_obs={
            frozenset({'dNch_deta', 'dET_deta', 'dN_dy'}),
        }
):
    """
    Estimate a covariance matrix for the given system and pair of observables,
    e.g.:

    >>> cov('PbPb2760', 'dN_dy', 'pion', 'dN_dy', 'pion')
    >>> cov('PbPb5020', 'dN_dy', 'pion', 'dNch_deta', None)

    For each dataset, stat and sys errors are used if available.  If only
    "summed" error is available, it is treated as sys error, and `stat_frac`
    sets the fractional stat error.

    Systematic errors are assumed to have a Gaussian correlation as a function
    of centrality percentage, with correlation length set by `sys_corr_length`.

    If obs{1,2} are the same but subobs{1,2} are different, the sys error
    correlation is reduced by `cross_factor`.

    If obs{1,2} are different and uncorrelated, the covariance is zero.  If
    they are correlated, the sys error correlation is reduced by
    `cross_factor`.  Two different obs are considered correlated if they are
    both a member of one of the groups in `corr_obs` (the groups must be
    set-like objects).  By default {Nch, ET, dN/dy} are considered correlated
    since they are all related to particle / energy production.

    """
    def unpack(obs, subobs):
        dset = data[system][obs][subobs]
        yerr = dset['yerr']

        try:
            stat = yerr['stat']
            sys = yerr['sys']
        except KeyError:
            stat = dset['y'] * stat_frac
            sys = yerr['sum']

        return dset['x'], stat, sys

    x1, stat1, sys1 = unpack(obs1, subobs1)
    x2, stat2, sys2 = unpack(obs2, subobs2)

    if obs1 == obs2:
        same_obs = (subobs1 == subobs2)
    else:
        # check if obs are both in a correlated group
        if any({obs1, obs2} <= c for c in corr_obs):
            same_obs = False
        else:
            return np.zeros((x1.size, x2.size))

    # compute the sys error covariance
    C = (
        np.exp(-.5*(np.subtract.outer(x1, x2)/sys_corr_length)**2) *
        np.outer(sys1, sys2)
    )

    if same_obs:
        # add stat error to diagonal
        C.flat[::C.shape[0]+1] += stat1**2
    else:
        # reduce correlation for different observables
        C *= cross_factor

    return C


def print_data(d, indent=0):
    """
    Pretty print the nested data dict.

    """
    prefix = indent * '  '
    for k in sorted(d):
        v = d[k]
        k = prefix + str(k)
        if isinstance(v, dict):
            print(k)
            print_data(v, indent + 1)
        else:
            if k.endswith('cent'):
                v = ' '.join(
                    str(tuple(int(j) if j.is_integer() else j for j in i))
                    for i in v
                )
            elif isinstance(v, np.ndarray):
                v = str(v).replace('\n', '')
            print(k, '=', v)


if __name__ == '__main__':
    print_data(data)
