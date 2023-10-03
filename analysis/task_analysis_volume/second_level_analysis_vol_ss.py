import sys
import os
import pandas as pd
import numpy as np
import glob
import json
import nibabel 
import nilearn
from nilearn.glm.contrasts import compute_fixed_effects
from nilearn.masking import intersect_masks
import matplotlib.pyplot as plt




#function to run second level fixed effects model to combine runs

def run_second_level_model(sub,task,ses):
    memory = './nilearn_cache' #change if desired
    space='MNI152NLin6Asym' #change if desired
    
    #grab masks needed from the two runs (could be run1 and run2, or run1 and run3 or run2 and run3)
    masks = glob.glob(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-*{space}*brain_mask.nii.gz')
    
    #make intersection of the two masks
    #threshold=1 corresponds to keeping the intersection of all masks, whereas threshold=0 is the union of all masks
    mask = intersect_masks(masks, threshold=1, connected=True)

        
    #check that both runs exist for a sub, task, ses
    #TO DO: technically should also check that the events.tsv exists for them for cleanest check
    if len(masks) == 1:
        print(f'Only one run exists for subject {sub}, task {task} and session {ses}. Second level analysis not possible.')
        return

    
    #find all contrasts for the two runs
    effect_size_list=glob.glob(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-*_contrast-*_effect_size.nii.gz')

    #create contrast to effect size path directory
    effect_size_dict={}
    for path in effect_size_list:
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
        effect_variance_list=glob.glob(f'../../../derivatives/task_analysis_volume/first_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_run-*_contrast-{contrast}*_effect_variance.nii.gz')
        
        #TO DO: figure out meaning of precision weighted
        #note that results are made up of these: fx_results = [fixed_fx_contrast, fixed_fx_variance, fixed_fx_tstat]
        fx_results = compute_fixed_effects(effect_size_list, effect_variance_list, mask, precision_weighted=False)

        output_dict[contrast]=fx_results 
    
    return output_dict




def save_second_level_outputs(sub,task,ses):
    print(f'Now processing subject {sub} for task {task} in session {ses}')
    output_dict = run_second_level_model(sub,task,ses)
    
    
    #check that output does not exist
    output = glob.glob(f'../../../derivatives/task_analysis_volume/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-*_effect_size_fx.nii.gz')
    if output: #checks if list is not empty
        print(f'At least partial 1st level output exists for subject {sub}, session {ses} and task {task}. This must be deleted for before generating new output.')
        return 
    
    
    for contrast in output_dict.keys():
        #calculate contrast maps for all contrasts
        contrast_output = output_dict[contrast]
        
        #create paths to output dir if not exist
        derivatives_path = '../../../derivatives'
        nilearn_output_path = os.path.join(derivatives_path, 'task_analysis_volume','second_level',f'sub-{sub}',f'ses-{ses}',f'task-{task}')
        if not os.path.isdir(nilearn_output_path):
            os.makedirs (nilearn_output_path)
            
        #save contrast maps to files
        #note the order of fixed effects outputs: fx_results = [fixed_fx_contrast, fixed_fx_variance, fixed_fx_tstat]
        contrast_output[0].to_filename(f'../../../derivatives/task_analysis_volume/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-{contrast}_effect_size_fx.nii.gz')
        contrast_output[1].to_filename(f'../../../derivatives/task_analysis_volume/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-{contrast}_effect_variance_fx.nii.gz')
        contrast_output[2].to_filename(f'../../../derivatives/task_analysis_volume/second_level/sub-{sub}/ses-{ses}/task-{task}/sub-{sub}_ses-{ses}_task-{task}_rec-unco_contrast-{contrast}_stat_fx.nii.gz')
    print(f'{sub} {task} {ses}: processing done')
        
    return




#loop through all sessions and runs that exist for a participant and task that are provided as input to this script   
def main(sub,task):
    space='MNI152NLin6Asym'
    ses_list=['baseline','1year']
    
    #loop through all potential sessions
    for ses in ses_list:
        
        #find all subjects who have task niftis for this task, session
        masks = glob.glob(f'../../../derivatives/ses-{ses}/sub-{sub}/ses-{ses}/func/sub-{sub}*task-{task}*rec-unco*run-*{space}*brain_mask.nii.gz')
        
        if (len(masks)>0):
            print('in main getting started')
            save_second_level_outputs(sub,task,ses)
            print(f'DONE: finished 2nd level GLM of {sub} for {task} during session {ses}')
        
        else:
            print(f'WARNING: could not find data for 2nd level GLM of {sub} for {task} during session {ses}')
                
    print('finishing up main')




if __name__ == "__main__":
    try:
        sub = sys.argv[1]
        task = sys.argv[2]
    except IndexError:
        print("You did not specify both a subject and a task (subject is first argument, task is second)")
        sys.exit(1)
    main(sub,task)


