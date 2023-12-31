#!/bin/tcsh
#PBS -N pcamal 
#PBS -l nodes=10:ppn=2,walltime=1:00:00 

# clean out module defaults (for now anyway)
module purge

# load in the mpi module needed.
module load mpi/mpich2-1.0.5_p4_gcc-3.4.6

#for mpich use:
module load mpiexec
set mpiexec_ops="-pernode"

#for ompi their mpiexec.
#mpiexec_ops='-bynode'

# qsub 
module load torque
cd /home/pppebay/share-src/pcamal/examples/bisected_cyl_2_procs
mpiexec -transform-hostname='s/$/-ib/g' pcamal_proto 2 bisected_cyl.pcamal.g meshed_bisected_cyl

