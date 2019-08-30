#!/usr/bin/bash

job_id = $1

for design_pt in {0..99}; do
    echo "catting $design_pt"
    cd CAT-RESULTS
    mkdir $job_id
    cat $SCRATCH/results/$job_id/$design_pt/[0-9]*.dat >> CAT-RESULTS/$job_id/$design_pt.dat
done
