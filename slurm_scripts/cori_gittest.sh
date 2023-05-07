#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=00:30:00 
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=haswell

#SBATCH -A m1641
#SBATCH --ntasks-per-node=1
#SBATCH -J github5_IC
#SBATCH -o output/github/github5_IC.o
#SBATCH -e output/github/github5_IC.e
#SBATCH --mail-user=wade.cappa@wsu.edu

# # module use /global/common/software/m3169/perlmutter/modulefiles
# module use /global/common/software/m3169/cori/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

mpirun -n 5 ./build/release/tools/mpi-greedi-im -i /global/cfs/cdirs/m1641/network-data/Binaries/github_binary.txt -w -k 16 -p -d IC -e 0.13 -o output/github/github5_IC.json --run-streaming=true --reload-binary -u --alpha=0.25
