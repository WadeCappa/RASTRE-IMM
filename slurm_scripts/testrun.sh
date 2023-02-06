#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 12:00:00 
#SBATCH -q overrun
#SBATCH -N 1024
#SBATCH --ntasks-per-node=1
#SBATCH -J Orkut1024
#SBATCH -o output/orkut/Orkut1024.o
#SBATCH -e output/orkut/Orkut1024.e
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
module load openmpi
#module load cudatoolkit/11.0

mpirun -n 1024 ./build/release/tools/mpi-greedi-im -i test-data/orkut_small.txt -w -k 100 -p -d IC -e 0.13 -o Orkut1024_results.json
