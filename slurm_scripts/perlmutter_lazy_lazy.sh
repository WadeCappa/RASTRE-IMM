#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -N 128
#SBATCH --ntasks-per-node=1
#SBATCH -J m128_lazy_wikipedia
#SBATCH -o /global/cfs/cdirs/m1641/network-results/algo-comparisons/wikipedia/m128_lazy_wikipedia.o
#SBATCH -e /global/cfs/cdirs/m1641/network-results/algo-comparisons/wikipedia/m128_lazy_wikipedia.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wade.cappa@wsu.edu

# module use /global/common/software/m3169/perlmutter/modulefiles
module use /global/common/software/m3169/perlmutter/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#module load PrgEnv-cray
# module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.24.3
module load cray-mpich
module load cray-libsci
#module load openmpi
#module load cudatoolkit/11.0

srun -n 128 ./build/release/tools/mpi-greedi-im -i /global/cfs/cdirs/m1641/network-data/Binaries/wikipedia_binary.txt -w -k 100 -p -d IC -e 0.13 -o /global/cfs/cdirs/m1641/network-results/algo-comparisons/wikipedia/m128_lazy_wikipedia.json --run-streaming=false --reload-binary -u 