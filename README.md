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

2. cifti smoothing using hcp tools 

3. first, second (for mid and sst tasks given they have two runs) and group level modeling for nifti (volumetric) data

4. first, second (for mid and sst tasks given they have two runs) and group level modeling for cifti (grayordinate) data

5. fsleyes visualization for nifti/volumetric results

6. final visualizations using fsleyes outputs and pillow for nifti/volumetric results

7. final visualizations using nilearn and pillow for cifti/grayordinate results


