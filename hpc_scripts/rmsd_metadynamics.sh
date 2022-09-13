#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 168:00

# job name (default = name of script file)
#BSUB -J "rmsd_metadynamics"
# Set output file
#BSUB -o  %J.out

# Set error file
#BSUB -e %J.stderr

# Specify node group
#BSUB -q gpuqueue

# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=96]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"


source ~/.bashrc
conda activate enhanced-sampling-tmem175
python ../metadynamics/rmsd_force.py -i ../systems/system00 -r ../systems/system01 -o ../../enhanced-sampling-tmem175-data/rmsd_metadynamics -y ../metadynamics/tic1_0.7.yaml
#python ../metadynamics/rmsd_force.py -i ../systems/system00 -r ../systems/system01 -o ../../enhanced-sampling-tmem175-data/rmsd_metadynamics -y ../metadynamics/cas.yaml
