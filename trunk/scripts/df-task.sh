#!/usr/bin/bash

# Sleep depending on task ID - TASK_ID*5 seconds
sleep $((TASK_ID * 5))
date 
echo Starting task $TASK_ID

# Invoke the NVT Df script:
python /home/andrewj/Projects/LatticeSwitch/scripts/nvt-df-run.py 0.58 2

