#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 168:00
#
# Set output file
#BSUB -o  md.26.0.out
#
# Set error file
#BSUB -eo md.26.0.stderr
#
# Specify node group
#BSUB -m "lu-gpu lv-gpu ld-gpu lt-gpu lp-gpu lg-gpu boson"

#
# nodes: number of nodes 
#BSUB -n 1 -R "rusage[mem=96]"

#
# job name (default = name of script file)
#BSUB -J "26 0"
source ~/.bashrc
conda activate simenv2
python []
