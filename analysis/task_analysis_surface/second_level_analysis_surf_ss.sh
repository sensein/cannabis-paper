#!/bin/bash

#SBATCH -p gablab
#SBATCH --time=01:00:00
#SBATCH --mem=20GB
#SBATCH -c 14
#SBATCH -n 1
#SBATCH -J second_level_analysis
#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL
#SBATCH -x node[041,054-060,100-115]

# grab these from submission script
task=$1
args=($@)
subjs=(${args[@]:1}) #:1 needed because argument 1 is the task, which we need to exclude


# index slurm array to grab subject
subject=${subjs[${SLURM_ARRAY_TASK_ID}]}
subject=${subject:4:5}

echo Submitted subject: ${subject}

#activate python env needed for the script
source activate /om2/user/dclb/.miniconda/envs/nilearn_older_v

# default command
cmd="python second_level_analysis_surf_ss.py $subject $task"

printf "Command:\n${cmd}\n"

# run it
eval $cmd
