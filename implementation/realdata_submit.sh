#!/bin/sh

fold=$1

R_script="real_cancer"

module load intel/2018.2 openmpi/2.1.3 r/3.5.0

R CMD BATCH "--args fold=$fold " $SCRATCH/Microbiome_classification/updated_code/${R_script}.R $SCRATCH/Microbiome_classification/updated_code/out/${R_script}_$fold.Rout

##submit
#for fold in {1..10}; do sbatch -N 1 --ntasks-per-node=40 -t 59:00 -o $SCRATCH/Microbiome_classification/updated_code/out/fold${fold}.out -J fold${fold} $SCRATCH/Microbiome_classification/updated_code/realdata_submit.sh $fold; done;
