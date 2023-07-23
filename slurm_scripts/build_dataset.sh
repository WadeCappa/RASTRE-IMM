#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 24:00:00
#SBATCH -q preempt
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J building_friendster_LT
#SBATCH -o /global/homes/w/wadecap/building_friendster_LT.o
#SBATCH -e /global/homes/w/wadecap/building_friendster_LT.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wade.cappa@wsu.edu

# module use /global/common/software/m3169/perlmutter/modulefiles
module use /global/common/software/m3169/perlmutter/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=128
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

srun -n 1 ~/ripples/build/release/tools/dump-graph -i /global/cfs/cdirs/m1641/network-data/test_data/com-friendster.ungraph.txt --distribution uniform -d LT -o /global/homes/w/wadecap/friendster_LT_binary.txt --scale-factor 0.1 --dump-binary
