#!/bin/bash
#SBATCH -J test_event            # Job name
#SBATCH -p skx-normal      # Queue (partition) name
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 05:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=everett.165@osu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A TG-PHY180035    # Allocation name (STARTUP ALLOC)
#SBATCH -o slurm/out-%j
#SBATCH -e slurm/err-%j


module load intel/18.0.2 cmake/3.7.1 gsl boost hdf5 eigen impi python3

source ../prepare_compilation_stampede2_2.sh
#source ../../jf_prepare.sh

module list

inputdir=../input-config
job=$SLURM_JOB_ID
#ntasks=$(( SLURM_JOB_NUM_NODES * SLURM_CPUS_ON_NODE ))
ntasks=1

srun run-events --nevents 1 --rankvar SLURM_PROCID --rankfmt "{:0${#ntasks}d}" --logfile $SCRATCH/$job.log --tmpdir=$SCRATCH/ --startdir=$inputdir $job.dat
