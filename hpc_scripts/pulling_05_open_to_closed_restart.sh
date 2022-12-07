#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 48:00

# Set output file
#BSUB -o  pulling_05_open_to_closed_restart.out

#BSUB -J "open_to_closed"

# Set error file
#BSUB -e pulling_05_open_to_closed_restart.stderr

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
-i ../systems/system00 \
-e ../systems/eBDims/open_to_closed_superposed_to_final.pdb \
-o ../../enhanced-sampling-tmem175-data/pulling/05.open_to_closed \
-p ../sim_params/pulling02.yaml \
-y ../pulling_methods/cas.yaml \
-r ../../enhanced-sampling-tmem175-data/pulling/05.open_to_closed/frame40/ \
--restart_from_frame 40
