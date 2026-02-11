#!/bin/bash

N_STARS=${1:-100}
R_SCALE=${2:-100.0}
RHO_0=${3:-0.001}

echo "Plotting Rotation Curve:"
echo "Number of Stars   : $N_STARS"
echo "Scale Radius      : $R_SCALE"
echo "DM Density        : $RHO_0"

python src/plot_rotation_curve.py --N_STARS $N_STARS --R_SCALE $R_SCALE --RHO_0 $RHO_0
