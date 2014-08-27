#!/bin/bash
#
#PBS -N RealSpeAna
#PBS -j oe
#PBS -o RealSpecAnaJob.log
#PBS -V
#PBS -l nodes=1:ppn=1,walltime=02:00:00

cd $PBS_O_WORKDIR

#echo "Start Analysis of COSY data"
./RealSpecAna &> RealSpecAnaResults.log

