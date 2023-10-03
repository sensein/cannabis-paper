#!/bin/bash

#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL
#SBATCH -p normal

#SBATCH -J second_level_analysis_submission_script

# provide the task to be analyzed after sbatch script_name.sh
base=../../..
task=$1

# go to dicom directory and grab all subjects we want to convert (in this case all that have session baseline)
pushd ${base} > /dev/null
subjs=($(ls sub-* -d | sed -r 's|/[^/]+$||'))
popd > /dev/null

# take the length of the array
# this will be useful for indexing later
len=$(expr ${#subjs[@]} - 1)


echo Spawning ${#subjs[@]} sub-job

# submit subject to fmriprep  processing
sbatch --array=0-$len second_level_analysis_vol_ss.sh $task ${subjs[@]}
