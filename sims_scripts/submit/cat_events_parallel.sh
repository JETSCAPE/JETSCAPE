#!/usr/bin/bash 

#the job_id is the first command line argument 

cd catted_results

mkdir $1
cd $1

mkdir Events
cd Events

mkdir main
cd ..

cd ..

cd ..

for design_pt in {0..499}; do
    echo "catting $design_pt"
    cat $SCRATCH/results/$1/$design_pt/[0-9]*.dat >> catted_results/test_parallel/Events/main/$design_pt.dat &
done

#wait til all processes have finished
wait

echo "Finished catting results! Goodbye"

#zip -r catted_results/$1.zip $1
