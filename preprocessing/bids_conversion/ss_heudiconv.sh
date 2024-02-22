#!/bin/bash
#SBATCH -p gablab
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH -J heudiconv
#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL

# grab these from submission script
input=$1
ses=$2
args=($@)
subjs=(${args[@]:2})

# set output for conversion
outdir=../../..

# save location of code directory
codedir=.

# index slurm array to grab subject
subject=${subjs[${SLURM_ARRAY_TASK_ID}]}

echo Submitted subject: ${subject}

module add openmind/singularity/3.6.3
SING_IMG=/om2/user/dclb/containers/imaging/heudiconv_0.11.3.sif

# default command
cmd="singularity run -B ${input}:/in -B ${outdir}:/out -B ${codedir}:/code ${SING_IMG} -d /in/{subject}/${ses}/MR* -o /out -f /code/heuristic.py -c dcm2niix -s ${subject} -ss ${ses} -b --minmeta -g accession_number"

# dry run command
#cmd="singularity run -B ${input}:/in m -B ${outdir}:/out ${SING_IMG} -d in/{subject}/${ses}/MR* -o /out -f convertall -c none -s ${subject} -ss ${ses} --minmeta -g accession_number"

printf "Command:\n${cmd}\n"

# run it
eval $cmd
