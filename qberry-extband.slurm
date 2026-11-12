#!/bin/bash
#SBATCH --partition=A1,A2,A3,A4,B1,B2,B3,B4                # select partition (A1, A2, A3, A4, B1, B2, or B3) 
#SBATCH --time=3:00:00               # set time limit in HH:MM:SS 
#SBATCH --nodes=1                     # number of nodes 
#SBATCH --ntasks-per-node=8          # number of processes per node (for MPI) 
#SBATCH --cpus-per-task=1             # OMP_NUM_THREADS (for openMP) 
#SBATCH --job-name=qberry                # job name #SBATCH --output="error.%x"
#SBATCH --qos=short                   # short for jobs up to 3 days,medium 7 days,long 14 days
                                      # standard output and error are redirected to 

                                      # <job name>_<job ID>.out # for OpenMP jobs export
                                      # if srun -n not working, try mpirun -np
#OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
set -e
OMP_NUM_THREADS= 1 
MKL_NUM_THREADS= 1

main=`pwd`

# v2qgt executable
module load oneapi/compiler/2022.0.2
module load oneapi/mpi/2021.5.1
module load oneapi/mkl/2022.0.2
module load oneapi/oclfpga/2022.0.1
v2qgt=/home/ytqian/soft/v2qgt/v2qgt_v5/v2qgt.x

# Target Kramers pair
NN=597
MM=598

# Active window: 109-112
BMIN=597
BMAX=599

# Derivative directions for plane A: Gamma-T and W-T
# G = (0,0,0)
# T = (0,0,0.5)
# W = (0.25,0.25,0.25)
#
# d1 = G - T = (0,0,-0.5)
# d2 = W - T = (0.25,0.25,-0.25)
export QGT_DVEC1="1.0 0.0  0.0"
export QGT_DVEC2="0.0 1.0  0.0"

export QGT_EXT_BMIN=$BMIN
export QGT_EXT_BMAX=$BMAX

for i in {42..100}
do
  file="file$i"

  if [ -d "$file" ]
  then
     cd "$file" || exit 1
     echo "$file"

     # ==================================================
     # all external bands
     # ==================================================
     export QGT_EXT_MODE=all

     mpirun -np $SLURM_NTASKS $v2qgt -nn $NN -mm $MM -kubo 7 -qm 0 -kp 1 -s 2 -abelian 2 -Qnm 1 > berryout_all

     # ==================================================
     # near active window only
     # For NN=109, MM=110, BMIN=109, BMAX=112,
     # near means external bands l = 111,112
     # ==================================================
     export QGT_EXT_MODE=near

     mpirun -np $SLURM_NTASKS $v2qgt -nn $NN -mm $MM -kubo 7 -qm 0 -kp 1 -s 2 -abelian 2 -Qnm 1 > berryout_near

     # ==================================================
     # remote bands only
     # remote means bands outside 109-112
     # ==================================================
     export QGT_EXT_MODE=remote

     mpirun -np $SLURM_NTASKS $v2qgt -nn $NN -mm $MM -kubo 7 -qm 0 -kp 1 -s 2 -abelian 2 -Qnm 1 > berryout_remote

     sleep 2

     cd "$main" || exit 1
  fi
done
