import os


def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where

    allowed template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """

    rest_moco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-rest_rec-moco_run-{item:1d}_bold')
    rest_unco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-rest_rec-unco_run-{item:1d}_bold')
    spin_task_AP=create_key('{bids_subject_session_dir}/fmap/sub-{subject}_{session}_acq-task_dir-AP_run-{item:1d}_epi')
    spin_task_PA=create_key('{bids_subject_session_dir}/fmap/sub-{subject}_{session}_acq-task_dir-PA_run-{item:1d}_epi')
    spin_diff=create_key('{bids_subject_session_dir}/fmap/sub-{subject}_{session}_acq-diff_dir-AP_run-{item:1d}_epi')
    dwi=create_key('{bids_subject_session_dir}/dwi/sub-{subject}_{session}_acq-{acq}_dwi')
    t1=create_key('{bids_subject_session_dir}/anat/sub-{subject}_{session}_run-{item:1d}_T1w')
    t2=create_key('{bids_subject_session_dir}/anat/sub-{subject}_{session}_run-{item:1d}_T2w')
    nback_moco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-nback_rec-moco_run-{item:1d}_bold')
    nback_unco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-nback_rec-unco_run-{item:1d}_bold')
    sst_moco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-sst_rec-moco_run-{item:1d}_bold')
    sst_unco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-sst_rec-unco_run-{item:1d}_bold')
    mid_moco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-mid_rec-moco_run-{item:1d}_bold')
    mid_unco=create_key('{bids_subject_session_dir}/func/sub-{subject}_{session}_task-mid_rec-unco_run-{item:1d}_bold')

    info = {rest_moco: [], rest_unco: [], spin_task_PA: [], spin_task_AP: [], spin_diff: [], dwi: [], t1: [], t2: [], nback_moco: [], nback_unco: [], sst_moco: [], sst_unco: [], mid_moco: [], mid_unco: []}
    

    #data = create_key('run{item:03d}')
    #info = {data: []}
    #last_run = len(seqinfo)

    for s in seqinfo:
        """
        The namedtuple `s` contains the following fields:

        * total_files_till_now
        * example_dcm_file
        * series_id
        * dcm_dir_name
        * unspecified2
        * unspecified3
        * dim1
        * dim2
        * dim3
        * dim4
        * TR
        * TE
        * protocol_name
        * is_motion_corrected
        * is_derived
        * patient_id
        * study_description
        * referring_physician_name
        * series_description
        * image_type
        """

        x,y,sl,nt = (s.dim1, s.dim2, s.dim3, s.dim4)
        

        if (sl == 256) and (nt == 1) and ('T1_p2' in s.protocol_name):
            info[t1].append(s.series_id)
        elif (sl == 256) and (nt == 1) and ('T2_SPACE' in s.protocol_name):
            info[t2].append(s.series_id)
        elif (nt == 74) and ('ep2d_diff_sms' in s.protocol_name):  
            if ('DFC' in s.image_type):
                info[dwi].append({'item': s.series_id, 'acq': 'DFC'})
            else:
                info[dwi].append({'item': s.series_id, 'acq': 'ORIG'})
        elif (nt == 1) and ('_topup_' in s.protocol_name):
            if ('task' in s.protocol_name):
                if ('ap' in s.protocol_name):
                    info[spin_task_AP].append(s.series_id)
                elif ('pa' in s.protocol_name):
                    info[spin_task_PA].append(s.series_id)
            elif ('diff' in s.protocol_name):
                info[spin_diff].append(s.series_id)
        elif (sl == 69) and ('_resting' in s.protocol_name) and (nt == 243):
            if s.is_motion_corrected:
                info[rest_moco].append(s.series_id)
            elif not s.is_motion_corrected:
                info[rest_unco].append(s.series_id)
        elif (sl == 69) and ('_MID' in s.protocol_name) and (nt == 215):
            if s.is_motion_corrected:
                info[mid_moco].append(s.series_id)
            elif not s.is_motion_corrected:
                info[mid_unco].append(s.series_id)   
        elif (sl == 69) and ('_SST' in s.protocol_name) and (nt == 257):
            if s.is_motion_corrected:
                info[sst_moco].append(s.series_id)
            elif not s.is_motion_corrected:
                info[sst_unco].append(s.series_id) 
        elif (sl == 69) and ('_Nback' in s.protocol_name) and (nt == 278):
            if s.is_motion_corrected:
                info[nback_moco].append(s.series_id)
            elif not s.is_motion_corrected:
                info[nback_unco].append(s.series_id)
        else:
            pass

    for key, items in info.items():
        print(items)

    return info
