import os
import sys
import glob
import json
import re

import pandas as pd
import numpy as np

import nibabel as nib
import nilearn

from nilearn.glm.first_level import FirstLevelModel
from nilearn.plotting import plot_design_matrix
from nilearn.plotting import plot_stat_map, plot_anat, plot_img
from nilearn.maskers import NiftiMasker
from nilearn import image as nimg

import matplotlib.pyplot as plt




#get confounds for any subject and task
def get_confounds(sub,task,ses,run):
    print('in get_confounds getting started')
    #load file with confounds timeseries
    all_confounds = pd.read_csv(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_desc-confounds_timeseries.tsv', sep = '\t')
    
    #rigid body motion
    motion_params = ['trans_x', 'trans_y', 'trans_z','rot_x','rot_y','rot_z']
    
    #individual col with single 1 for timepoint of motion
    motion_outliers = [col for col in all_confounds.columns if 'motion_outlier' in col]  
    
    #for low freq signal drift
    #cannot include this and high-pass temp filter bc already removes low freq fluc
    #required if using aCompCor (or tCompCor)
    cosine_regressors = [col for col in all_confounds.columns if 'cosine' in col] 
    
    #these can be adjusted to be from the combined wm csf, for example
    #doesn't make sense to use csf and wm signal regression if using these according to fmriprep documentation
    #6 is rule of thumb; can pick diff number or specific amount of variance explained
    num_a_comp_cors=6
    a_comp_cors = []
    for i in range(num_a_comp_cors):
        a_comp_cors.append('a_comp_cor_{:02d}'.format(i))
    
    #for now leave out crown regressors but can be added using the following code later on    
    #same for edge comp cors
    num_edge_comp_cors=6
    edge_comp_cors = []
#     for i in range(num_edge_comp_cors):
#         edge_comp_cors.append('edge_comp_{:02d}'.format(i))
    
#     #OR instead get all edge comp cors
#     edge_comp_cors = [col for col in all_confounds.columns if 'edge_comp' in col]
        
    #we need to filter out non-steady state volumes if using cosine regressors, ICA AROMA and CompCor regressors...    
    non_steady_state_regressors = [col for col in all_confounds.columns if 'non_steady_state' in col]
    
    #merge desired confounds
    selected_confounds = all_confounds[['framewise_displacement']+motion_params+motion_outliers+cosine_regressors+a_comp_cors+edge_comp_cors+non_steady_state_regressors].copy()

    #get rid of nas in first row of derivative and framewise displacement cols
    for col in selected_confounds.columns:
        if ('derivative' in col) or ('framewise_displacement' in col):
            if pd.isna(selected_confounds[col][0]):
                selected_confounds[col][0]=0
    print('nuisance regressors used:')
    print(selected_confounds.columns)
    print('finishing up get_confounds')
    return selected_confounds




#run first level glm model
def run_first_level_model(sub,task,ses,run):
    print('in run_first_level_model getting started')
    space='MNI152NLin6Asym' #change if desired
      
    #grab nifti, events and mask files
    nifti = glob.glob(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-{run}*{space}*preproc*nii.gz')[0]
    events = glob.glob(f'../../../sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-{run}*events.tsv')[0]
    mask = glob.glob(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-{run}*{space}*brain_mask.nii.gz')[0]
    
    events = pd.read_csv(events, sep = '\t')
    
    #events file needs to be modified to reflect combined image+arrow stimuli for the sst task
    if task == 'sst':
        events=events.iloc[::2].reset_index(drop=True)
        events['duration']=events['duration'] + 1.5
        events['trial_type'] = events['trial_type'].replace({'N_Go_image': 'Go', 
                                                             'MJ_Go_image': 'Go', 
                                                             'N_UnsuccStop_image': 'UnsuccStop', 
                                                             'MJ_UnsuccStop_image': 'UnsuccStop',
                                                             'N_SuccStop_image': 'SuccStop', 
                                                             'MJ_SuccStop_image': 'SuccStop'})
        
    #call function to get confounds
    selected_confounds=get_confounds(sub,task,ses,run)
    
    #get TR
    task_json = open(f"../../../task-{task}_bold.json")
    task_json=json.load(task_json)
    TR=task_json['RepetitionTime']
    
    
    #smoothe and mask nifti
    #smoothing fwhm selected to be 4 since it is about twice the voxel resolution
    nifti_masker = NiftiMasker(
        mask_img = mask,
        standardize=False,
        standardize_confounds=False,
        mask_strategy="epi",
        smoothing_fwhm=4,
        t_r = TR)
    
    nifti_masked_smoothed_2d = nifti_masker.fit_transform(nifti)
    nifti_masked_smoothed = nifti_masker.inverse_transform(nifti_masked_smoothed_2d)
    
    
    #get median grayordinate val
    HCP_smoothed_cifti = glob.glob(f'../../../derivatives/HCP_smoothing/sub-{sub}/ses-{ses}/smoothed_sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_space-fsLR_den-91k_bold.dtseries.nii')[0]
    HCP_smoothed_cifti_nimg = nimg.load_img(HCP_smoothed_cifti)
    HCP_smoothed_cifti_signal_unscaled = HCP_smoothed_cifti_nimg.get_fdata(dtype='f4')
    median_cifti_across_time_space = np.median(HCP_smoothed_cifti_signal_unscaled)
    
    #scale nifti timeseries
    masked_smoothed_nifti_nimg = nimg.load_img(nifti_masked_smoothed)
    masked_smoothed_nifti_signal_unscaled = masked_smoothed_nifti_nimg.get_fdata(dtype='f4')
    min_across_time_space = np.min(nifti_masked_smoothed_2d) #taking min of masker objects avoids getting 0 as min due to mask
    masked_smoothed_scaled_nifti_signal = 1000*(masked_smoothed_nifti_signal_unscaled-min_across_time_space)/(median_cifti_across_time_space-min_across_time_space)

    #get the affine matrix
    nifti_image = nib.load(nifti)
    affine_matrix = nifti_image.affine

    #create nifti image from the masked, smoothed and scaled data array
    masked_smoothed_scaled_nifti_nimg = nib.Nifti1Image(masked_smoothed_scaled_nifti_signal, affine=affine_matrix)

    #choose middle slice as reference since fmriprep default for slice timing correction
    #already accounting for drift and high pass filter with cosine regressors, so not needed
    #spm + derivative + disperson should account for undershoot and slight variances in HRF across ppl/regions
    #smoothing was already done above, so set to None here
    #mask has to be reapplied since otherwise the area around the brain would now be negative due to having subtracted the min from the 0-masked voxels
    glm = FirstLevelModel(t_r = TR, mask_img=mask, \
        slice_time_ref=0.5, smoothing_fwhm=None, hrf_model='spm + derivative + dispersion', drift_model=None, \
        high_pass=None, drift_order=4, standardize=False, signal_scaling=False, noise_model='ar1', \
        minimize_memory=True, verbose=0, n_jobs=-2)

    #fit glm
    fitted_glm = glm.fit(masked_smoothed_scaled_nifti_nimg, events=events, confounds=selected_confounds)
   
    print('finishing up run_first_level_model')
    
    return [fitted_glm, selected_confounds]




#create all desired contrasts
def make_localizer_contrasts(design_matrix,task,selected_confounds):
    print('in make_localizer_contrasts getting started')

    # first generate canonical contrasts based on the design matrix columns 
    contrast_matrix = np.eye(design_matrix.shape[1])
    canonical_contrasts = dict([(column, contrast_matrix[i])
                      for i, column in enumerate(design_matrix.columns)])
    
    #initialize dictionary of contrasts desired for analysis of respective task
    final_contrasts={}
    
    #initialize list of complex contrasts, which are combinations of the canonical contrasts
    task_contrasts = []
    
    
    #complex contrasts for nback task (must include a '-' or '+' and be made up of canonical contrasts)
    if task == 'nback':
        task_contrasts = ['twoback-zeroback']
    
    #complex contrasts for mid task (must include a '-' or '+' and be made up of canonical contrasts)
    if task == 'mid':
        task_contrasts = ['HiRewCue-NeuCue', #high reward anticipation -- paper A2, ABCD
                          'HiRewCue-LoRewCue', #high vs. low reward anticipation -- paper A3, ABCD
                          'HiRewCue+LoRewCue-NeuCue', #combined reward anticipation -- paper A1
                          'LoRewCue-NeuCue', #combined reward anticipation -- ABCD
                          
                          'HiLossCue-NeuCue', #high loss anticipation -- paper A5, ABCD
                          'HiLossCue-LoLossCue', #high vs. low loss anticipation -- ABCD
                          'HiLossCue+LoLossCue-NeuCue', #combined loss anticipation
                          'LoLossCue-NeuCue', #combined loss anticipation -- ABCD

                          'HiRewCue-HiLossCue', #high reward vs. high loss anticipation
                          'HiRewCue+LoRewCue-HiLossCue-LoLossCue', #combined reward vs. loss anticipation

                          
                          'HiWin-NeuHit', #high reward outcome cp. to neutral hit -- paper O6
                          'HiWin+LoWin-NeuHit', #combined reward outcome cp. to neutral hit
                          
                          'HiWin-HiNoWin', #high reward outcome cp. to high reward miss
                          'HiWin+LoWin-HiNoWin-LoNoWin', #combined reward outcome cp. to combined reward miss -- ABCD
                          
                          'HiLoss-NeuMiss', #high loss cp. to neutral miss
                          'HiLoss+LoLoss-NeuMiss', #combined loss cp. to neutral miss
                          
                          'HiLoss-AvoidHiLoss', #high loss cp. to high avoid loss
                          'HiLoss+LoLoss-AvoidHiLoss-AvoidLoLoss', #combined loss cp. to combined avoid loss -- ABCD
                          
                          'HiLoss-NeuHit', #high loss cp. to neutral hit -- paper O7
                          'HiLoss+LoLoss-NeuHit', #combined loss cp. to neutral hit
                          
                          'HiWin-HiLoss', #high reward outcome cp. to high loss
                          'HiWin+LoWin-HiLoss+LoLoss', #combined reward outcome cp. to combined loss
                         ]
    
    if task == 'sst':
        #image and arrow as combined regressors 
        #MJ and N as combined regressors
        task_contrasts = ['SuccStop-Go', #correct inhibition
                          'UnsuccStop-Go', #incorrect inhibition
                          'UnsuccStop-SuccStop' #successful inhibitory control  
                         ]
    
    
    #for loop that creates complex contrasts based on the string names from the task_contrasts list created above
    for task_contrast in task_contrasts:
        #split complex contrast string into a list of canonical contrasts separated by '-' or '+' 
        #'-' and '+' are maintained in the list due to the () around the deliminator specified in re.split
        events = re.split('([\-\+])', task_contrast)
        
        #check if all canonical contrasts needed for a given complex contrast were present in a subject's data
        #note that sometimes a participant eg. never had a miss on a type of trial, so the miss regressor doesn't exist 
        #can comment out print statements to see which contrasts are being generated vs. for which events are missing
        if len(set(events).intersection(set(canonical_contrasts.keys()))) < np.ceil(len(events)/2):
            continue
        
        #set up placeholder for array containing matrix for the complex contrast
        complex_contrast = 0
        #set up how the canonical contrasts will be added/subtracted to create the complex contrast
        math_sign = '+'
        plus_multiplier = 1
        minus_multiplier = 1
        #this only applies to the mid task, since the nback and sst tasks' complex contrasts are only made up of two canonical contrasts
        #if mid list consists of 5 elements (i.e. 3 canonical contrasts and 2 math symbols), there is always 1 '+' and 1 '-' (in addition to the implied '+ at the beginning)
        #in that case, we need to multiply the 2 canonical constrasts after the '+' and implied '+' by 0.5
        if len(events) == 5:
            plus_multiplier = 0.5
        #if mid list consists of 7 elements (i.e. 4 canonical contrasts and 3 math symbols), there is always 1 '+' and 2 '-' (in addition to the implied '+ at the beginning)
        #in that case, we need to multiply the 4 canonical constrasts by 0.5
        if len(events) == 7:
            minus_multiplier = 0.5
            plus_multiplier = 0.5
            
        #do the complex contrast calculation
        for el in events:
            #math_sign starts of with '+' as set up above since the first canonical contrast has an implied '+' at the beginning (added to the 0 that we used to initialize the complex contrast)
            if el in ['+','-']:
                math_sign = el
            #if current list element isn't a math sign, then it's a canonical contrast and we add/subtract based on the math sign we set
            #we also multiply the canonical contrast by 1 (default) or 0.5 (certain mid complex contrasts)
            else:
                if math_sign == '+':
                    complex_contrast = complex_contrast + plus_multiplier*canonical_contrasts[el]
                if math_sign == '-':
                    complex_contrast = complex_contrast - minus_multiplier*canonical_contrasts[el]
        #at the end we store the created complex contrast        
        final_contrasts[task_contrast] = complex_contrast       
                
    #add a nuisance regressor contrast to the final set of contrasts that includes all nuisance regressors
    #firs set up zeros matrix of correct shape
    final_contrasts['nuisance_regressors']=np.zeros(design_matrix.shape[1])
    #this iteratively adds a 1 in place of all the confounds in the zeros matrix
    for confound in selected_confounds.columns:
        final_contrasts['nuisance_regressors']=final_contrasts['nuisance_regressors']+canonical_contrasts[confound]
    
    #add big win anticipation minues implicit baseline as contrast for mid
    if task == 'mid':
        final_contrasts['HiRewCue'] = canonical_contrasts['HiRewCue'] #high reward anticipation vs. baseline -- paper A4
    
    #return dictionary of all complex contrasts
    print('finishing up make_localizer_contrasts')
    return final_contrasts




#store first level GLM results
def save_first_level_outputs(input_list):
    print('in save_first_level_outputs getting started')
    space='MNI152NLin6Asym'
    sub,task,ses,run = input_list
    print(f'Now processing subject {sub} for task {task} in session {ses} run {run}')
    
    #check that input exists
    nifti = glob.glob(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-{run}*{space}*preproc*nii.gz')
    if not nifti: #checks if list is empty
        print(f'1st level input does not exist for subject {sub}, session {ses}, run {run} and task {task}. Cannot generate output.')
        return
    
    #check that output does not exist
    output = glob.glob(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_contrast-*_effect_size.nii.gz')
    if output: #checks if list is not empty
        print(f'At least partial 1st level output exists for subject {sub}, session {ses}, run {run} and task {task}. This must be deleted for before generating new output.')
        return    
    
    #run first level model 
    fitted_glm, selected_confounds = run_first_level_model(sub,task,ses,run)         
    
    #define contrasts for task
    contrasts = make_localizer_contrasts(fitted_glm.design_matrices_[0],task,selected_confounds)

    for contrast in contrasts.keys():
        #calculate contrast maps for all contrasts
        #F test to see if any nuisance regressor in contrast is significant
        if contrast == 'nuisance_regressors':
            contrast_output = fitted_glm.compute_contrast(contrasts[contrast], output_type='all',stat_type='F')
        #t test to see if specified individual contrast is significant
        else:   
            contrast_output = fitted_glm.compute_contrast(contrasts[contrast], output_type='all',stat_type='t')
        
        #rename mid contrasts to simpler names
        #set up dictionary with simpler names
        replace_mid_contrasts_dict = {
                                      'HiRewCue+LoRewCue-NeuCue':'RewCue-NeuCue', #reward anticipation
                                      'HiLossCue+LoLossCue-NeuCue':'LossCue-NeuCue', #loss anticipation
                                      'HiRewCue+LoRewCue-HiLossCue-LoLossCue': 'RewCue-LossCue', #combined reward vs. loss anticipation
                                      'HiRewCue':'HiRewCue-Baseline', #high reward anticipation vs. baseline

                                      'HiWin+LoWin-NeuHit':'Win-NeuHit', #reward outcome cp. to neutral hit
                                      'HiWin+LoWin-HiNoWin-LoNoWin':'Win-NoWin', #reward outcome cp. to reward miss
                                      'HiLoss+LoLoss-NeuMiss':'Loss-NeuMiss', #loss cp. to neutral miss
                                      'HiLoss+LoLoss-AvoidHiLoss-AvoidLoLoss':'Loss-AvoidLoss', #loss cp. to high loss
                                      'HiLoss+LoLoss-NeuHit':'Loss-NeuHit', #combined loss cp. to neutral hit
                                      'HiWin+LoWin-HiLoss+LoLoss':'Win-Loss', #combined reward outcome cp. to combined loss
                                      }
        #replace name of contrast if current contrast found in dictionary as a key, otherwise keeps current name
        contrast = replace_mid_contrasts_dict.get(contrast, contrast) 
        
        #create paths to output dir if not exist
        derivatives_path = '../../../derivatives'
        nilearn_output_path = os.path.join(derivatives_path, 'task_analysis_volume','first_level',f'sub-{sub}',f'ses-{ses}',f'task-{task}')
        if not os.path.isdir(nilearn_output_path):
            os.makedirs (nilearn_output_path)
        
        #save contrast maps to files
        contrast_output['effect_size'].to_filename(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_contrast-{contrast}_effect_size.nii.gz')
        contrast_output['effect_variance'].to_filename(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_contrast-{contrast}_effect_variance.nii.gz')
        contrast_output['z_score'].to_filename(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_contrast-{contrast}_z_score.nii.gz')
        contrast_output['stat'].to_filename(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_contrast-{contrast}_stat.nii.gz')
        contrast_output['p_value'].to_filename(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-{run}_contrast-{contrast}_p_value.nii.gz')
    
    print('finishing up save_first_level_outputs')
    return




#loop through all sessions and runs that exist for a participant and task that are provided as input to this script   
def main(sub,task):
    space='MNI152NLin6Asym'
    ses_list=['baseline','1year']
    run_list=['1','2','3']
    
    #loop through all potential sessions and runs
    for ses in ses_list:
        for run in run_list:
            
            #check that sub has task nifti and events.tsv file for this task, session and run
            nifti = (glob.glob(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-{run}*{space}*preproc*nii.gz'))
            events = (glob.glob(f'../../../sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-{run}*events.tsv'))
            #checks that nifti and events files exist and that events file isn't empty (if empty but nifti exists then this nifti is from a run that was repeated and should not be analyzed)
            if (len(nifti)>0) and (len(events)>0) and pd.read_csv(events[0],sep='\t').empty == False:
                print('in main getting started')
                input_list = sub,task,ses,run
                save_first_level_outputs(input_list)
                print(f'DONE: finished 1st level GLM of {sub} for {task} during session {ses} run {run}')
            else:
                print(f'WARNING: could not find data for 1st level GLM of {sub} for {task} during session {ses} run {run}')
                
    print('finishing up main')
    

if __name__ == "__main__":
    try:
        sub = sys.argv[1]
        task = sys.argv[2]
    except IndexError:
        print("You did not specify both a subject and a task (subject is first argument, task is second)")
        sys.exit(1)
    main(sub,task)



