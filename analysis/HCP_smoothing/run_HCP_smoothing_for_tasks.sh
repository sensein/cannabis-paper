#!/bin/bash

#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL

#SBATCH -J HCP_submission_script


tasks=(nback sst mid rest)

runs=(1 2 3)

for task in ${tasks[@]}
do
    echo "Start HCP wb smoothing of all surfaces from task $task"
    for run in ${runs[@]}
    do
        echo "Run $run"
        sbatch ./HCP_smoothing_by_task.sh $task $run
    done
    echo "Finished HCP wb smoothing of all surfaces from task $task"
done

echo "HCP wb smoothing finished for all tasks and runs!"



