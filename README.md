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


### Preprocessing:

1. Heudiconv-based BIDS conversion

2. fmriprep 23.0.1 preprocessing pipeline

3. mriqc-based quality metric generation
  
4. Events file population with appropriate event timings and durations


### Analysis:

1. Demographics and cannabis characteristics by group table generation

* Run the following jupyter notebooks in order:
  * _analysis/behavioral_analysis/cannabis_table.ipynb_
  * _analysis/behavioral_analysis/demographics_table.ipynb_


2. cifti smoothing using hcp tools 

* Run the following using SLURM job scheduler:
  * _analysis/HCP_smoothing/run_HCP_smoothing_for_tasks.sh_
    * this will call for each task: _analysis/HCP_smoothing/HCP_smoothing_by_tasks.sh_


3. first, second (for mid and sst tasks given they have two runs) and group level modeling for nifti (volumetric) data

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


4. first, second (for mid and sst tasks given they have two runs) and group level modeling for cifti (grayordinate) data

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




