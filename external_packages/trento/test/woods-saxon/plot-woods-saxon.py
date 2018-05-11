#!/usr/bin/env python3

"""
Compile the C++ Woods-Saxon generator, run it, pipe the results into this
script, and plot.
"""

import os.path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


def woods_saxon(r, R, a):
    return r*r/(1. + np.exp((r-R)/a))


def main():
    basename = 'generate-woods-saxon'
    cxxname = '{}.cxx'.format(basename)
    if not os.path.exists(basename) or (
            os.path.getmtime(cxxname) > os.path.getmtime(basename)):
        print('compiling C++ Woods-Saxon generator')
        subprocess.check_call(
            ['g++', '-std=c++11', '-Wall', '-Wextra', '-O3', '-march=native',
             cxxname, '-o'+basename])

    R = 6.5
    a = 0.5
    N = int(1e7)

    print('generating Woods-Saxon numbers')
    with subprocess.Popen(['./'+basename, str(R), str(a), str(N)],
                          stdout=subprocess.PIPE) as proc:
        samples = np.fromiter(proc.stdout, dtype=float, count=N)

    print('plotting')
    plt.hist(samples, bins=200, histtype='step', normed=True)

    rmax = samples.max()
    r = np.linspace(0, rmax, 1000)
    norm = integrate.quad(woods_saxon, 0, rmax, args=(R, a))[0]
    plt.plot(r, woods_saxon(r, R, a)/norm)

    plt.xlim(0, rmax)
    plt.show()


if __name__ == "__main__":
    main()
