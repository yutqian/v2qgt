#!/bin/bash 
#SBATCH --partition=A1                # select partition (A1, A2, A3, A4, B1, B2, or B3) 
#SBATCH --time=6-24:00:00               # set time limit in HH:MM:SS 
#SBATCH --nodes=1                     # number of nodes 
#SBATCH --ntasks-per-node=20          # number of processes per node (for MPI) 
#SBATCH --cpus-per-task=1             # OMP_NUM_THREADS (for openMP) 
#SBATCH --job-name=mos                # job name #SBATCH --output="error.%x"
                                      # standard output and error are redirected to 
                                      # <job name>_<job ID>.out # for OpenMP jobs export
                                      # if srun -n not working, try mpirun -np
OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# load modules if needed

module load oneapi/mkl/2024.1
module load oneapi/mpi/2021.12
module load oneapi/oclfpga/2024.1.0
module load oneapi/compiler/2024.1.0

vasp=/ioss/ytqian/soft/vasp-wann1.2/bin/vasp_ncl
##mpirun -np $SLURM_NTASKS $vasp  > vasp.out


srun -n 16  $vasp  > vasp.out
