#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 1:30:00
#SBATCH -q regular
#SBATCH -N 256
#SBATCH --ntasks-per-node=1
#SBATCH -J orkut_mpi_greedi_im_256
#SBATCH -o output/orkut/Orkut256.o
#SBATCH -e output/orkut/Orkut256.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wade.cappa@wsu.edu

# module use /global/common/software/m3169/perlmutter/modulefiles
module use /global/common/software/m3169/perlmutter/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#module load PrgEnv-cray
# module load python/3.9-anaconda-2021.11
module load gcc/11.2.0
module load cmake/3.24.3
module load openmpi
#module load cudatoolkit/11.0

srun -n 256 --mpi=pmi2 ./build/release/tools/mpi-greedi-im -i test-data/orkut_small.txt -w -k 50 -p -d IC -e 0.13 -o output/orkut/orkut256_result.json