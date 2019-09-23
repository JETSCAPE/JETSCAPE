#!/usr/bin/bash 

cd catted_results

#this is the destination folder 
mkdir $1
cd $1

#the next five clis are the five job IDs of design points that all will be catted together

mkdir Events
cd Events

mkdir main
cd ..

cd ..

cd ..

#loop over all design points 
for design_pt in {0..499}; do
    echo "catting $design_pt"
    #run each command in background to fully utilize multiprocessing
    {
    cat $SCRATCH/results/$2/$design_pt/[0-9]*.dat >> catted_results/$1/Events/main/$design_pt.dat
    cat $SCRATCH/results/$3/$design_pt/[0-9]*.dat >> catted_results/$1/Events/main/$design_pt.dat
    cat $SCRATCH/results/$4/$design_pt/[0-9]*.dat >> catted_results/$1/Events/main/$design_pt.dat
    cat $SCRATCH/results/$5/$design_pt/[0-9]*.dat >> catted_results/$1/Events/main/$design_pt.dat
    cat $SCRATCH/results/$6/$design_pt/[0-9]*.dat >> catted_results/$1/Events/main/$design_pt.dat
    } &
done

#wait till all processes have finished to print success
wait

echo "All results have been concatenated! Goodbye"

#zip -r catted_results/$1.zip $1
