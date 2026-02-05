#!/bin/bash

#Default arguments: N_STARS, R_SCALE, RHO_0, DT, STEPS
N_STARS=${1:-100}
R_SCALE=${2:-100.0}
RHO_0=${3:-0.001}
DT=${4:-0.005} # Default time step increased significantly
STEPS=${5:-2} #Calculate 2 physics steps for every 1 frame drawn
DEBUG=${6:-0} # 0 = OFF, 1 = ON

echo "Launching Simulation with:"
echo "__________________________"
echo "Number of Stars   : $N_STARS"
echo "Scale Radius      : $R_SCALE"
echo "DM Density        : $RHO_0"
echo "Time Step (dt)    : $DT"
echo "Steps per Frame   : $STEPS"
echo "Debug Mode        : $DEBUG"
echo "__________________________"

CMD = "python tests/visualize.py \
    --N_STARS $N_STARS \
    --R_SCALE $R_SCALE \
    --RHO_0 $RHO_0 \
    --DT $DT \
    --STEPS $STEPS"

# check for debug flag
if [ "$DEBUG" -eq "1" ]; then
    CMD="$CMD --DEBUG_TREE"
fi

# run the command
$CMD
