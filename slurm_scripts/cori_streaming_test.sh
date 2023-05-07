#!/bin/bash

#SBATCH --qos=debug
#SBATCH --time=00:30:00 
#SBATCH --nodes=33
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=haswell

#SBATCH -A m1641
#SBATCH --ntasks-per-node=1
#SBATCH -J Orkut33_streaming
#SBATCH -o output/orkut/Orkut33_test.o
#SBATCH -e output/orkut/Orkut33_test.e
#SBATCH --mail-user=wade.cappa@wsu.edu

# # module use /global/common/software/m3169/perlmutter/modulefiles
# module use /global/common/software/m3169/cori/modulefiles

#OpenMP settings:
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

mpirun -n 33 ./build/release/tools/mpi-greedi-im -i /global/cfs/cdirs/m1641/network-data/Binaries/orkut_small_binary.txt -w -k 100 -p -d IC -e 0.13 -o Orkut33_test.json --run-streaming=true --epsilon-2=0.1 --reload-binary
