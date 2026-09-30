#!/bin/bash
#PBS -r n
#PBS -m n
#PBS -N aiida-None
#PBS -V
#PBS -o _scheduler-stdout.txt
#PBS -e _scheduler-stderr.txt
#PBS -q MONARIS
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
cd "$PBS_O_WORKDIR"


export QE_MPI_RANKS=88

'/home/trevizam/programs/qe-7.6-environ-3.1.1-gcc14-openmpi5-mkl2023-full/aiida-wrappers/pw-wrapper.x' '-nk' '8' '-in' 'aiida.in'  > 'aiida.out'
