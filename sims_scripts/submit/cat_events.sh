#!/usr/bin/bash

for i in {0..99}; do
    echo "catting $i" 
    cat $SCRATCH/results/*/$i/[0-9]* > CAT-RESULTS/$i.dat
done
