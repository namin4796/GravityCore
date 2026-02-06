#!/bin/bash


CONFIG=${1:-"default_config.json"}
DEBUG=${2:-1} # 0 = OFF, 1 = ON

echo "__________________________"
echo "Launching Simulation with:"
echo "__________________________"


CMD="python tests/visualize.py --config $CONFIG"

# check for debug flag
if [ "$DEBUG" -eq "1" ]; then
    CMD="$CMD --DEBUG_TREE"
fi

# run the command
$CMD
