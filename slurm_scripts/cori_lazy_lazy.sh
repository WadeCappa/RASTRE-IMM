#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=00:30:00 
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=haswell

#SBATCH -A m1641
#SBATCH --ntasks-per-node=1
#SBATCH -J Github4_lazy_lazy
#SBATCH -o output/github/Github4_lazy_lazy.o
#SBATCH -e output/github/Github4_lazy_lazy.e
#SBATCH --mail-user=wade.cappa@wsu.edu

# # module use /global/common/software/m3169/perlmutter/modulefiles
# module use /global/common/software/m3169/cori/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

mpirun -n 4 ./build/release/tools/mpi-greedi-im -i test-data/githubSmall.txt -w -k 100 -p -d IC -e 0.13 -o Github4_lazy_lazy.json --run-streaming=false