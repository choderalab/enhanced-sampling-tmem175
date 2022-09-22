#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 24:00

# Set output file
#BSUB -o  %J.out

#BUSB -J "pulling"

# Set error file
#BSUB -e %J.stderr

# Specify node group
#BSUB -q gpuqueue
#BSUB -R select[gpu_model0!='GeForceGTX1080']

# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"


source ~/.bashrc
conda activate enhanced-sampling-tmem175
python ../pulling_methods/pulling_sequential.py -i ../systems/system00 -e ../systems/eBDims/open_to_closed_superposed_to_final.pdb -o ../../enhanced-sampling-tmem175-data/pulling/open_to_closed -p ../sim_params/pulling.yaml -y cas.yaml
