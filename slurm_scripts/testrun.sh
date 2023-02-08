#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 12:00:00 
#SBATCH -q overrun
#SBATCH -N 64
#SBATCH --ntasks-per-node=1
#SBATCH -J Orkut64
#SBATCH -o /global/homes/w/wadecap/ripples/output/orkut/Orkut64.o
#SBATCH -e /global/homes/w/wadecap/ripples/output/orkut/Orkut64.e
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

mpirun -n 64 ./build/release/tools/mpi-greedi-im -i test-data/orkut_small.txt -w -k 100 -p -d IC -e 0.13 -o Orkut64_results.json
