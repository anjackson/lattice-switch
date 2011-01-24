#!/usr/bin/bash
#$ -A d04-ajac
#$ -l h_rt=71:59:
### -l h_rt=47:00:
### -l h_rt=00:29:
#$ -pe hpc 8
### -t 1-4:1
#$ -cwd
#$ -S /usr/bin/bash
# Export shell variables:
#$ -V

# We want Grid Engine to send mail
# when the job begins
# and when it ends.

#$ -M anj@anjackson.net
#$ -m e

# We want to name the file for the standard output
# and standard error.

#$ -o df-para.out -j y

mprun -np 8 /home/andrewj/Projects/LatticeSwitch/scripts/mpi_task_launcher /home/andrewj/Projects/LatticeSwitch/scripts/df-task.sh

