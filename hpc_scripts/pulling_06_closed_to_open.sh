#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 60:00

# Set output file
#BSUB -o  %J.out

#BSUB -J "closed_to_open"

# Set error file
#BSUB -e %J.stderr

# Specify node group
#BSUB -q gpuqueue
#BSUB -R select[gpu_model0!='NVIDIAGeForceGTX1080']
#BSUB -R select[gpu_model0!='NVIDIAGeForceGTX1080Ti']

# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=8]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"

source ~/.bashrc
conda activate enhanced-sampling-tmem175
python ../pulling_methods/pulling_sequential.py \
-i ../systems/system01 \
-e ../systems/eBDims/closed_to_open_superposed_to_final.pdb \
-o ../../enhanced-sampling-tmem175-data/pulling/06.closed_to_open \
-p ../sim_params/pulling02.yaml \
-y ../pulling_methods/cas.yaml

#python ../pulling_methods/pulling_sequential.py \
#-i ../systems/system01 \
#-e ../systems/eBDims/closed_to_open_superposed_to_final.pdb \
#-o ../../enhanced-sampling-tmem175-data/pulling/04.closed_to_open_restart \
#-p ../sim_params/pulling02.yaml \
#-y ../pulling_methods/cas.yaml \
#-r ../../enhanced-sampling-tmem175-data/pulling/04.closed_to_open/frame40/ \
#--restart_from_frame 40
