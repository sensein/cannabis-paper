#!/bin/bash

#SBATCH -p gablab
#SBATCH --time=96:00:00
#SBATCH --mem=35GB
#SBATCH --cpus-per-task=16
#SBATCH -J fp_baseline
#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL
#SBATCH -x node[054-060,100-115]

# grab these from submission script
bids_dir=$1
session=$2
args=($@)
subjs=(${args[@]:2}) #:2 needed because argument 1 is the input directory and argument 2 is the session, which we need to exclude

# set output for conversion
outdir=${bids_dir}/derivatives/ses-${session}
export SINGULARITYENV_TEMPLATEFLOW_HOME=/home/dclb/.cache/templateflow

# index slurm array to grab subject
subject=${subjs[${SLURM_ARRAY_TASK_ID}]}

echo Submitted subject: ${subject}

#set up singularity container
module add openmind/singularity/3.6.3
SING_IMG=/om2/user/dclb/containers/imaging/fmriprep_23.0.1.sif

# assign working directory
scratch=/om2/scratch/Sun/$(whoami)/fmriprep_work/${session}/${subject}
mkdir -p ${scratch}

#set fs license
export SINGULARITYENV_FS_LICENSE=/om2/user/dclb/imaging/fs_license_linux.txt

#set up
unset $PYTHONPATH
set -e

# default command
cmd="singularity run -e -B /om2 -B ${scratch}:/workdir -B ${bids_dir}:/input -B ${outdir}:/out ${SING_IMG} /input /out participant --participant_label ${subject} -w /workdir --skip-bids-validation --use-aroma --cifti-output 91k -v -v --mem_mb 40000 --output-spaces MNI152NLin6Asym:res-2 anat --force-bbr --bids-filter-file /input/code/preprocessing/fmriprep_MM_${session}/${session}_unco_filter.json --nprocs 16 --omp-nthreads 8 --bids-database-dir /out/bids_db/"

printf "Command:\n${cmd}\n"

#remove an IsRunning files from freesurfer
rm -f ${outdir}/sourcedata/freesurfer/${subject}/scripts/IsRunning*

# run it
eval $cmd
