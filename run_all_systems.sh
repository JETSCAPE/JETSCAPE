#!/bin/bash

p3 make-ee-bins.py -d /data/rjfgroup/rjf01/cameron.parker/tunes/$1/parameters.txt -n $1-tune/LEP -p
p3 make-pp-bins.py -d /data/rjfgroup/rjf01/cameron.parker/tunes/$1/parameters.txt -n $1-tune/LHC2760 -p
p3 make-pp-bins.py -d /data/rjfgroup/rjf01/cameron.parker/tunes/$1/parameters.txt -n $1-tune/LHC13000 -p
p3 make-pp-bins.py -d /data/rjfgroup/rjf01/cameron.parker/tunes/$1/parameters.txt -n $1-tune/RHIC -p