#!/bin/bash

#SBATCH -p gablab
#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL

#SBATCH -J heudiconv_submission_script

input=../../../sourcedata/dicom/

# go to dicom directory and grab all subjects we want to convert
pushd ${input} > /dev/null
subjs=($(ls *_* -d)) 
#subjs=($(ls MM_122 -d)) # can instead select just one for debugging here, for example MM_122
popd > /dev/null

# take the length of the array
# this will be useful for indexing later
len=$(expr ${#subjs[@]} - 1)


ses=(baseline 1year)

echo Spawning ${#subjs[@]} sub-jobs for each session.

# submit subject to heudiconv processing
for s in ${ses[@]}; do
	sbatch --array=0-$len ss_heudiconv.sh $input $s ${subjs[@]}
done
