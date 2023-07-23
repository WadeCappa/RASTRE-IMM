#!/bin/bash

#SBATCH -A m1641
#SBATCH -C cpu
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -N 512
#SBATCH --ntasks-per-node=1
#SBATCH -J diimm_512_IC_friendster
#SBATCH -o /global/homes/w/wadecap/results/jobs/friendster/diimm_512_IC_friendster.o
#SBATCH -e /global/homes/w/wadecap/results/jobs/friendster/diimm_512_IC_friendster.e
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

srun -n 512 ./build/release/tools/mpi-imm -i /global/cfs/cdirs/m1641/network-data/Binaries/friendster_IC_binary.txt -w -k 100 -p -d IC -e 0.13 -o /global/homes/w/wadecap/results/jobs/friendster/diimm_512_IC_friendster.json --reload-binary --DIiMM=true -u
