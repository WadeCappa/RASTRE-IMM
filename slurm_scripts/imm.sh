#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 00:30:00
#SBATCH -q regular
#SBATCH -N 256
#SBATCH --ntasks-per-node=1
#SBATCH -J imm_256_LT_livejournal
#SBATCH -o /global/homes/w/wadecap/strong_scaling/livejournal/imm_256_LT_livejournal.o
#SBATCH -e /global/homes/w/wadecap/strong_scaling/livejournal/imm_256_LT_livejournal.e
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


srun -n 256 ./build/release/tools/mpi-imm -i /global/cfs/cdirs/m1641/network-data/Binaries/livejournal_LT_binary.txt -w -k 100 -p -d LT -e 0.13 -o /global/homes/w/wadecap/strong_scaling/livejournal/imm_256_LT_livejournal.json --reload-binary 
