#!/usr/bin/bash

#for i in {0..200}; do
#    echo "catting $i" 
#    cat $SCRATCH/results/*/$i/[0-9]* > CAT-RESULTS/$i.dat
#done

for i in {0..49}; do
    echo "catting $i" 
    cat $SCRATCH/results-validation/4114650/$i/[0-9]* > CAT-RESULTS-VALID/$i.dat
done
