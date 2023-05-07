#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=00:30:00 
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=haswell

#SBATCH -A m1641
#SBATCH --ntasks-per-node=1
#SBATCH -J Orkut16_lazy_lazy
#SBATCH -o output/orkut/Orkut16_lazy_lazy.o
#SBATCH -e output/orkut/Orkut16_lazy_lazy.e
#SBATCH --mail-user=wade.cappa@wsu.edu

# # module use /global/common/software/m3169/perlmutter/modulefiles
# module use /global/common/software/m3169/cori/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

mpirun -n 16 ./build/release/tools/mpi-greedi-im -i test-data/orkut_small.txt -w -k 16 -p -d IC -e 0.13 -o Orkut16_lazy_lazy.json --run-streaming=false