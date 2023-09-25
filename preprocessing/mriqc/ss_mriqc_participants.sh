#!/bin/bash

#SBATCH -p normal
#SBATCH --time=3-23
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=16
#SBATCH -J mriqc_participants
#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL

# grab these from submission script
bids_dir=$1
args=($@)
subjs=(${args[@]:1}) #:1 needed because argument 1 is the input directory, which we need to exclude

myscratch=/om2/scratch/Sun/$(whoami)

# set output for conversion
outdir=${bids_dir}/derivatives/mriqc

# index slurm array to grab subject
subject=${subjs[${SLURM_ARRAY_TASK_ID}]}

echo Submitted subject: ${subject}

#set up singularity container
module add openmind/singularity/3.6.3
SING_IMG=/om2/user/dclb/containers/imaging/mriqc_0.16.1.sif

# assign working directory
scratch=${myscratch}/mriqc_work/${subject}
mkdir -p ${scratch}

#set up
unset $PYTHONPATH
set -e

# default command
cmd="singularity run -e -B /om2 -B ${scratch}:/workdir -B ${bids_dir}:/input:ro -B ${outdir}:/out ${SING_IMG} /input /out participant --participant_label ${subject} -w /workdir --nprocs 16 --omp-nthreads 8"

printf "Command:\n${cmd}\n"


# run it
eval $cmd
