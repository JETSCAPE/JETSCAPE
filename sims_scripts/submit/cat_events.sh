#!/usr/bin/bash

for i in {0..99}; do
    echo "catting $i" 
    cat $SCRATCH/results/$1/$i/[0-9]*.dat >> CAT-RESULTS/$i.dat
done
