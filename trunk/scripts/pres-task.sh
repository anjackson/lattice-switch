#!/usr/bin/bash

# Sleep depending on task ID - TASK_ID*5 seconds
sleep $((TASK_ID * 5))
date 
echo Starting task $TASK_ID

# Invoke the EOS script:
/home/andrewj/Projects/LatticeSwitch/scripts/pres-dens-run.py 

