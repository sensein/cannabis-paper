# No changes in functional brain activation during cognitive processes detected after 12 months of using cannabis for medical symptoms

Debbie Burdinski, Alisha Kodibagkar, Kevin Potter, Randi Schuster, A. Eden Evins, Jodi Gilman*, Satrajit Ghosh*

*co-supervised


## Code Execution Guidelines

### Directory structure:

Main directory: 
* Contains:
  * BIDS-organized nifti and cifti data with matching events files
  * _code_ directory
  * _derivatives_ directory

Code directory: 
* Structured like this GitHub repository
* Code is to be executed inside the folder in which a given script is located (to use relative paths)
* Contains:
  * _analysis_ directory with code for behavioral nback analysis, HCP-based smoothing, first, second and group linear modeling and demographics/cannabis characteristics table generation
  * _environment_setup_ directory with yml file for the needed python environment 
  * _preprocessing_ directory with code for heudiconv, fmriprep, events file population and mriqc

Derivatives directory:
* Saves outputs from preprocessing and analyses and the final visualizations
  * Contains:
    * _behavioral_ directory with results from the nback behavioral analysis
    * _demograpgics_ directory with demographics and cannabis characteristics by group tables
    * _HCP_smoothing_ directory with smoothed cifti files using the HCP pipeline's tool
    * _mriqc_ directory with mriqc metrics by scan
    * _ses-1year_ directory with fmriprep outputs from the one-year timepoint scans
    * _ses-baseline_ directory with fmriprep outputs from the baseline timepoint scans
    * _task_analysis_surface_ directory with first and second level outputs from the cifti-based analysis and intermediate and final visualizations
    * _task_analysis_volume_ directory with first and second level outputs from the nifti-based analysis and intermediate and final visualizations
   

### Code Requirements:

1. Create a conda environment using _environment_setup/cannabis_paper_env.yml_ and load it before running any python scripts
2. Download the following singularity containers (eg. using Dockerhub):
  * heudiconv: heudiconv_0.11.3.sif
  * fmriprep: fmriprep_23.0.1.sif
    * in order to run fmriprep successfully, note that you also need a valid fs_license.txt (change path accordingly)
  * mriqc: mriqc_0.16.1.sif
3. Install the following:
  * singularity: Singularity 3.9.5
  * hcp tools: Connectome Workbench v1.2.3

_Note that the path the containers and installations will need to be changed depending on where you keep yours!_

### Preprocessing:

1. Heudiconv-based BIDS conversion

* Run the following using SLURM job scheduler:
  * _preprocessing/bids_conversion/submit_job_array.sh_
    * this will call for each participant: _preprocessing/bids_conversion/ss_heudiconv.sh_
    * which will use: _preprocessing/bids_conversion/heuristic.py_
   
* Run the following python script to assign field maps to scans:
  * _preprocessing/bids_conversion/add_intendedfor_excep.py_


2. fmriprep preprocessing pipeline

* HC baseline: Run the following using SLURM job scheduler:
  * _preprocessing/fmriprep_HC_baseline/submit_job_array_baseline.sh_
    * this will call for each task: _preprocessing/fmriprep_HC_baseline/ss_fmriprep_baseline.sh_
    * which will use: _preprocessing/fmriprep_HC_baseline/baseline_unco_filter.json_
   
* MMC baseline: Run the following using SLURM job scheduler:
  * _preprocessing/fmriprep_MM_baseline/submit_job_array_baseline.sh_
    * this will call for each task: _preprocessing/fmriprep_MM_baseline/ss_fmriprep_baseline.sh_
    * which will use: _preprocessing/fmriprep_MM_baseline/baseline_unco_filter.json_
      
* MMC one year: Run the following using SLURM job scheduler:
  * _preprocessing/fmriprep_MM_1year/submit_job_array_1year.sh_
    * this will call for each task: _preprocessing/fmriprep_MM_1year/ss_fmriprep_1year.sh_
    * which will use: _preprocessing/fmriprep_MM_1year/1year_unco_filter.json_


3. mriqc-based quality metric generation

* Run the following using SLURM job scheduler:
  * _preprocessing/mriqc/submit_job_array_mriqc_participants.sh_
    * this will call for each task: _preprocessing/mriqc/ss_mriqc_participants.sh_

* Run the following jupyter notebooks in order:
  * _preprocessing/mriqc/create_mriqc_tsv.ipynb_


4. Events file population with appropriate event timings and durations

* Run the following jupyter notebooks in order:
  * _preprocessing/events_file_population/create_SST_events_files.ipynb_
  * _preprocessing/events_file_population/create_nback_events_files.ipynb_
  * _preprocessing/events_file_population/create_MID_events_files.ipynb_
    * note that an alternative script for MID to model both cue and feedback exists here: _preprocessing/events_file_population/create_MID_events_files_cue_and_feedback.ipynb_
    

### Analysis:

1. Demographics and cannabis characteristics by group table generation

* Run the following jupyter notebooks in order:
  * _analysis/behavioral_analysis/cannabis_table.ipynb_
  * _analysis/behavioral_analysis/demographics_table.ipynb_


2. Cifti smoothing using hcp tools 

* Run the following using SLURM job scheduler:
  * _analysis/HCP_smoothing/run_HCP_smoothing_for_tasks.sh_
    * this will call for each task: _analysis/HCP_smoothing/HCP_smoothing_by_tasks.sh_


3. First, second (for mid and sst tasks given they have two runs) and group level modeling for nifti (volumetric) data

* First level: Run the following using SLURM job scheduler:
  * Replace {task} with the task that you want to run the first level model for
  * _analysis/task_analysis_volume/submit_job_array_first_level.sh {task}_
    * this will call for each participant: _analysis/task_analysis_volume/first_level_analysis_ss.sh_
    * which will call for each participant: _analysis/task_analysis_volume/first_level_analysis_ss.py_

* Second level (only for mid and sst tasks given they have two runs): Run the following using SLURM job scheduler:
  * Replace {task} with the task that you want to run the first level model for
  * _analysis/task_analysis_volume/submit_job_array_second_level_vol.sh {task}_
    * this will call for each participant: _analysis/task_analysis_volume/second_level_analysis_vol_ss.sh_
    * which will call for each participant: _analysis/task_analysis_volume/second_level_analysis_vol_ss.py_

* Group level: Run the following jupyter notebook in order:
  * Set _task=_ in the block that loops through all relevant participants to the one you want to analyze
  * Note that you can use commenting to make modeling decisions, eg. outlier selection, sub-group analysis and covariates 
  * _analysis/task_analysis_volume/group_level_analysis_vol.ipynb_
    * this will save the effect size and stat maps per group/session/task at _../../../derivatives/task_analysis_volume/group_level/group-{group}/ses-{ses}/task-{task}_ to be used by the fsleyes visualization


4. First, second (for mid and sst tasks given they have two runs) and group level modeling for cifti (grayordinate) data

* First level: Run the following using SLURM job scheduler:
  * Replace {task} with the task that you want to run the first level model for
  * _analysis/task_analysis_surface/submit_job_array_first_level.sh {task}_
    * this will call for each participant: _analysis/task_analysis_surface/first_level_analysis_ss.sh_
    * which will call for each participant: _analysis/task_analysis_surface/first_level_analysis_ss.py_

* Second level (only for mid and sst tasks given they have two runs): Run the following using SLURM job scheduler:
  * Replace {task} with the task that you want to run the first level model for
  * _analysis/task_analysis_surface/submit_job_array_second_level_surf.sh {task}_
    * this will call for each participant: _analysis/task_analysis_surface/second_level_analysis_surf_ss.sh_
    * which will call for each participant: _analysis/task_analysis_surface/second_level_analysis_surf_ss.py_

* Group level: Run the following jupyter notebook in order:
  * Set _task=_ in the block that loops through all relevant participants to the one you want to analyze
  * Note that you can use commenting to make modeling decisions, eg. outlier selection, sub-group analysis and covariates 
  * _analysis/task_analysis_surface/group_level_analysis_surf.ipynb_
      * this will save the left/right hemispheres and flat maps as well as the coronal display per group/session/task at _../../../derivatives/task_analysis_surface/visualization/raw_figures_ to be used by the final pillow visualization


5. fsleyes visualization for nifti/volumetric results

* Run the following jupyter notebook in order:
  * Set _task=_ in the block that loops through all relevant participants to the one you want to analyze
  * _analysis/fsleyes_vis/fsleyes_visualization.ipynb_
    * this will save the fsleyes visualizations per group/session/task at _../../../derivatives/task_analysis_volume/visualization/raw_figures_ to be used by the final pillow visualization


6. final visualizations using fsleyes outputs and pillow for nifti/volumetric results

* Run the following jupyter notebook in order:
  * Set _task=_ in the block that loops through all relevant participants to the one you want to analyze
  * _analysis/task_analysis_volume/final_volume_plots.ipynb_
    * this will save the final visualizations per task at _../../../derivatives/task_analysis_volume/visualization/final_nifti_figures_ 


7. final visualizations using nilearn and pillow for cifti/grayordinate results

* Run the following jupyter notebook in order:
  * Set _task=_ in the block that loops through all relevant participants to the one you want to analyze
  * _analysis/task_analysis_surface/final_volume_plots.ipynb_
    * this will save the intermediate and final visualizations per task at _../../../derivatives/task_analysis_surface/visualization/intermediate_cifti_figures_ 
    * this will save the final visualizations per task at _../../../derivatives/task_analysis_surface/visualization/final_cifti_figures_ 




