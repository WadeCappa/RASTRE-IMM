#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=00:30:00 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=haswell

#SBATCH -A m1641
#SBATCH --ntasks-per-node=1
#SBATCH -J Github1
#SBATCH -o output/github/Github1.o
#SBATCH -e output/github/Github1.e
#SBATCH --mail-user=wade.cappa@wsu.edu

# # module use /global/common/software/m3169/perlmutter/modulefiles
# module use /global/common/software/m3169/cori/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

mpirun -n 1 ./build/release/tools/mpi-greedi-im -i test-data/githubSmall.txt -w -k 16 -p -d IC -e 0.13 -o Github1.json --run-streaming=true
