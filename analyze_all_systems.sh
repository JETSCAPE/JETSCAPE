#!/bin/bash

if [ -d "$1/LEP/" ]; then
    echo "LEP analysis"
    p3 run-ee-analysis.py $1/LEP/ -p
fi

if [ -d "$1/LHC13000/" ]; then
    echo "LHC 13000 analysis"
    p3 run-pp-analysis.py $1/LHC13000/ -p
fi

if [ -d "$1/LHC2760/" ]; then
    echo "LHC 2760 analysis"
    p3 run-pp-analysis.py $1/LHC2760/ -p
fi

if [ -d "$1/RHIC/" ]; then
    echo "RHIC analysis"
    p3 run-pp-analysis.py $1/RHIC/ -p
fi