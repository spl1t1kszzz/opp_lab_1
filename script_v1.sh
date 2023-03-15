#!/bin/sh

#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=10000m,place=scatter
#PBS -m n
#PBS -o out.txt
#PBS -e err.txt

cd $PBS_O_WORKDIR
echo
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
./parallel_v1.out
