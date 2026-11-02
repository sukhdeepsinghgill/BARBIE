#!/bin/bash
#SBATCH -J zbs            
#SBATCH -p medium
#SBATCH -o %j-out
#SBATCH -e %j-err
#SBATCH -N 2
#SBATCH --ntasks-per-node=40                                  
#SBATCH --cpus-per-task=1         
#SBATCH -t 50:00:00                

module load compiler/intel/2022/mpi/oneapi-2021.5.1 compiler/intel/2022/oneapi-2022.0.2
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/20ph91f04/local_install/lib/


#module load compiler/openmpi/4.0.2

#exe=./zbs_v3_mpi 384 137
exe=./zbs_fast_mpi 384 101

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#run the application
mpirun -n $SLURM_NTASKS $exe    # specify the application command-line options, if any, after $exe
