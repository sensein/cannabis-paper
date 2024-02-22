#!/bin/bash

#SBATCH -p gablab
#SBATCH --time=2:00:00
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=2
#SBATCH -J HCP_smoothing_by_task
#SBATCH -o log/%x-%A-%a.out
#SBATCH -e log/%x-%A-%a.err
#SBATCH --mail-user=dclb@mit.edu
#SBATCH --mail-type=ALL
#SBATCH -x node[054-060,100-115]



task=$1 #read in task to be smoothed
run=$2 #read in run of the task to be smoothed

WORKPATH="../../../derivatives"

OUTPATH="../../../derivatives/HCP_smoothing"
mkdir -p $OUTPATH
OUTPATH=`realpath $OUTPATH`

cd $WORKPATH

for ses in * ;
do
    [[ $ses != ses-* ]] && continue
    echo "smoothing cifti files from $ses"
    cd $ses
    for sub in * ;
    do
        input_cifti_starred_path="./$sub/$ses/func/*task-$task*run-$run*dtseries.nii"
        [[ $sub == *.html ]] && continue
        if [ ! -f $input_cifti_starred_path ]; then
            echo "file $input_cifti_starred_path not found!"
            continue
        fi
        echo "mkdir -p $OUTPATH/$sub/$ses"
        mkdir -p "$OUTPATH/$sub/$ses"
        input_cifti_path=`realpath $input_cifti_starred_path`
        input_cifti=(${input_cifti_path//// })
        input_cifti=${input_cifti[12]}
        output_cifti=smoothed_$input_cifti
        output_cifti_path=$OUTPATH/$sub/$ses/$output_cifti
        echo "running wb_command -cifti-smoothing"
        echo "input $input_cifti_path"
        echo "output $output_cifti_path"
        now=$(date)
        echo "$now"
        #command inputs in order: cifti file; surface kernel size in terms of sigma, volume kernel size in terms of sigma; direction to smooth along: COLUMN is used for dtseries files, while ROW would be used for dconn files; cifti output file; left surface file; right surface file
        wb_command -cifti-smoothing \
        $input_cifti_path \
        1.6986 1.6986 \
        COLUMN \
        $output_cifti_path \
        -left-surface ../../code/analysis/templates/tpl-fsLR_hemi-L_den-32k_sphere.surf.gii \
        -right-surface ../../code/analysis/templates/tpl-fsLR_hemi-R_den-32k_sphere.surf.gii
    done
    cd ..
done

