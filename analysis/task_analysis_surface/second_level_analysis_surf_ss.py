import sys
import os
import pandas as pd
import numpy as np
import glob
import json

import nilearn
import nibabel as nb

from nilearn.glm.contrasts import _compute_fixed_effects_params
from nilearn import image as nimg


#function to run second level fixed effects model to combine runs

def run_second_level_model(sub,task,ses):
    memory = 'nilearn_cache' #change if desired
    
    #find all contrasts for the two runs
    effect_size_list_all = glob.glob(f'../../../derivatives/task_analysis_surface/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-*_contrast-*_effect_size.dscalar.nii')
    
    #create contrast to effect size path directory
    effect_size_dict = {}
    for path in effect_size_list_all:
        contrast = path.split('contrast-')[-1].split('_effect')[0]
        if contrast not in effect_size_dict.keys():
            effect_size_dict[contrast]=[]
        effect_size_dict[contrast].append(path)
        
    #do fixed effects analysis of all available contrasts 
    output_dict={}
    for contrast in effect_size_dict.keys():
        effect_size_list = effect_size_dict[contrast]
        #check if it was possible to calculate a given contrast for both runs
        #if not, skip this contrast
        if len(effect_size_list) < 2:
            print(f'Only one {contrast} contrast for {sub} during {task} in {ses}. Second level analysis not possible.')
            continue
        effect_variance_list=glob.glob(f'../../../derivatives/task_analysis_surface/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-*_contrast-*_variance.dscalar.nii')
        
        #TO DO: figure out meaning of precision weighted
        #note that results are made up of these: fx_results = [fixed_fx_contrast, fixed_fx_variance, fixed_fx_tstat]
        fx_results = _compute_fixed_effects_params(
            np.squeeze(
                [nb.load(fname).get_fdata(dtype='f4') for fname in effect_size_list]
            ),
            np.squeeze(
                [nb.load(fname).get_fdata(dtype='f4') for fname in effect_variance_list]
            ),
            precision_weighted=False)
        
        output_dict[contrast]=fx_results 
    
    return output_dict



def save_second_level_outputs(sub,task,ses):
    print(f'Now processing subject {sub} for task {task} in session {ses}')
    output_dict = run_second_level_model(sub,task,ses)
    
    
    #check that output does not exist
    output = glob.glob(f'../../../derivatives/task_analysis_surface/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-*_effect_size_fx.dscalar.nii')
    if output: #checks if list is not empty
        print(f'At least partial 2nd level output exists for subject {sub}, session {ses} and task {task}. This must be deleted for before generating new output.')
        return 
    
    orig_dscalar_list = glob.glob(f'../../../derivatives/task_analysis_surface/first_level/sub-*/ses-{ses}/task-{task}/sub-*_ses-{ses}_task-{task}_rec-unco_run-*_contrast-*_effect_size.dscalar.nii')
    orig_orig_dscalar_nimg = nimg.load_img(orig_dscalar_list[0])
    time_axis, brain_model_axis = [orig_orig_dscalar_nimg.header.get_axis(i) for i in range(orig_orig_dscalar_nimg.ndim)]

    for contrast in output_dict.keys():
        #calculate contrast maps for all contrasts
        contrast_output = output_dict[contrast]
        dscalar_output_dict = create_dscalar(brain_model_axis, contrast_output) 
        
        #create paths to output dir if not exist
        derivatives_path = '../../../derivatives'
        nilearn_output_path = os.path.join(derivatives_path, 'task_analysis_surface','second_level',f'sub-{sub}',f'ses-{ses}',f'task-{task}')
        if not os.path.isdir(nilearn_output_path):
            os.makedirs (nilearn_output_path)
            
        #save contrast maps to files
        dscalar_output_dict['effect_size_fx'].to_filename(f'../../../derivatives/task_analysis_surface/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-{contrast}_effect_size_fx.dscalar.nii')
        dscalar_output_dict['effect_variance_fx'].to_filename(f'../../../derivatives/task_analysis_surface/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-{contrast}_effect_variance_fx.dscalar.nii')
        dscalar_output_dict['stat_fx'].to_filename(f'../../../derivatives/task_analysis_surface/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-{contrast}_stat_fx.dscalar.nii')
    
    print(f'{sub} {task} {ses}: processing done')
        
        
        
def create_dscalar(brain_model_axis, contrast_output):
    
    #note the order of fixed effects outputs: fx_results = [fixed_fx_contrast, fixed_fx_variance, fixed_fx_tstat]
    contrast_output_dict = {'effect_size_fx':contrast_output[0], 'effect_variance_fx':contrast_output[1], 'stat_fx':contrast_output[2]}
    
    dscalar_output_dict = {}
    
    for output_type,output_value in contrast_output_dict.items():
        #create new header
        scalar_axis = nb.cifti2.ScalarAxis([output_type])
        new_header = nb.Cifti2Header.from_axes([scalar_axis, brain_model_axis])
        
        #ensure dimensions fit
        output_value = np.atleast_2d(output_value)
    
        #create dscalar
        dscaler = nb.Cifti2Image(output_value, new_header)
        
        #save in dict
        dscalar_output_dict[output_type]=dscaler
    
    return dscalar_output_dict

    

#loop through all sessions and runs that exist for a participant and task that are provided as input to this script   
def main(sub,task):
    ses_list=['baseline','1year']
    
    #loop through all potential sessions
    for ses in ses_list:
        
        #find all subjects who have task niftis for this task, session
        first_level_outputs = glob.glob(f'../../../derivatives/task_analysis_surface/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-*_contrast-*_effect_size.dscalar.nii')
        
        if (len(first_level_outputs)>0):
            print('in main getting started')
            save_second_level_outputs(sub,task,ses)
            print(f'DONE: finished 2nd level GLM of {sub} for {task} during session {ses}')
        
        else:
            print(f'WARNING: could not find data for 2nd level GLM of {sub} for {task} during session {ses}')
                
    print('finishing up main')
    return



if __name__ == "__main__":
    try:
        sub = sys.argv[1]
        task = sys.argv[2]
    except IndexError:
        print("You did not specify both a subject and a task (subject is first argument, task is second)")
        sys.exit(1)
    main(sub,task)