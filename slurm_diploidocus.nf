#!/bin/bash -l

#SBATCH --job-name=diploidocus
#SBATCH --account=$PAWSEY_PROJECT
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --export=ALL

module load nextflow/23.10.0
module load minimap2/2.26--he4a0461_1
module load singularity/3.11.4-nompi
nextflow run main.nf --projectDir . --assembly path/to/assembly/file 
