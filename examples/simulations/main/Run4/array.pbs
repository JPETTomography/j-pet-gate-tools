#!/bin/sh

#PBS -q i3d
#PBS -l nodes=1:ppn=1
#PBS -N 192str_Run4
#PBS -V # All env. variables are passed into the cluster

cd ${PBS_O_WORKDIR}    # Take me to the directory where I launched qsub

Gate ./.Gate/main/main${PBS_ARRAYID}.mac

exit 0;
