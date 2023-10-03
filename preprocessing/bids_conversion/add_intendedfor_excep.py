#!/usr/bin/env python

# Adopted from Gregory Ciccarelli's orginal script of the same name

from bids import BIDSLayout
from glob import glob
import os
import os.path as op
import json
import logging
import sys


if op.exists('intendedfor.log'):
    raise RuntimeError('Already generated - delete log to rerun')

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='intendedfor.log',
                    filemode='w')

ROOT = '../../..'
layout = BIDSLayout(ROOT)

# key: fieldmap run ID
# values: corresponding functional scans

def process_fieldmaps(fmaps, subj, session):
    fmap_count = 0
    dup_runs = False
    desired_runs = ['1','1'] #runs to use from the AP task and PA task scans for the given session (always 1 for diff in this dataset) 
    if subj in ['MM286','MM334','HC012','MM122','MM240']: #note that even though HC_033 has dup runs in ses-baseline, the desired runs are both 1, and thus the subj is not included here
        desired_runs=set_runs(subj,session)
    for fmap in fmaps:
        targets=[]
        run = fmap.split('_run-')[-1][0:1]
        if fmap.split('_acq-')[-1][:4] == 'diff':
            fmap_count += 1
            targets = grab_targets(['dwi'], subj, session, 'diff')
        elif fmap.split('_acq-')[-1][:4] == 'task':
            if (run == desired_runs[0] and fmap.split('_dir-')[-1][:2] == 'AP'):
                fmap_count += 1
                targets = grab_targets(layout.get_tasks(), subj, session, 'task')
            elif (run == desired_runs[1] and fmap.split('_dir-')[-1][:2] == 'PA'):
                fmap_count += 1
                targets = grab_targets(layout.get_tasks(), subj, session, 'task')
            else:
                dup_runs = True
        else:
            logging.warning('{} session {} has incorrect fieldmap labeling. "intended for" was not added to the mislabeled fieldmap run. Fix and rerun'.format(subj,session))
        if targets:
            intended = {'IntendedFor': targets}
            add_metadata(fmap, intended)
    if dup_runs:
        logging.warning('{} session {} has duplicate fieldmaps. "intended for" was not populated for scans other than the indicated set of 3 fieldmaps.'.format(subj,session))
    if fmap_count < 3:
        logging.warning('{} session {} does not have all 3 types of fieldmap runs (AP/PA task and AP diff).'.format(subj,session))
    else:
        logging.info('{} session {} had "intended for" successfully added for the indicated set of 3 fieldmaps'.format(subj,session))

def set_runs(subj,session):
    desired_runs=['1','1']
    if (subj == 'MM286' and session == '1year') or (subj == 'HC012' and session == 'baseline'):
        desired_runs=['2','2']
    elif (subj == 'MM334' and session == '1year') or (subj == 'MM122' and session == 'baseline'):
        desired_runs = ['3','3']
    elif (subj == 'MM334' and session == 'baseline'):
        desired_runs=['3','2']
    elif (subj == 'MM240' and session == '1year'):
        desired_runs=['2','1']
    return desired_runs

def grab_targets(vals, subj, session, target):
    targets = []
    for val in vals:
        if target == 'task':
            targets += glob(op.join(ROOT, 'sub-{}'.format(subj), 'ses-{}'.format(session), 'func', '*{}*nii.gz'.format(val)))
        elif target == 'diff':
            targets += glob(op.join(ROOT, 'sub-{}'.format(subj), 'ses-{}'.format(session), 'dwi', '*{}*nii.gz'.format(val)))
    if not targets:
        logging.warning('{} session {} does not have any {} scans'.format(subj,session,target))
        return

    return relpath(targets, subj)

def relpath(targets, subj):
    splitter = op.join(ROOT, 'sub-{}/'.format(subj))
    return [x.split(splitter)[-1] for x in targets]

def load_json(filename):
    """ easy load of json dict """
    with open(filename, 'r') as fp:
        data = json.load(fp)
    return data

def add_metadata(fl, data, ind=4):
    """Adds dict items to exisiting json
    Parameters
    ----------
    fl : File (path to json file)
    data : dict (items to add)
    ind: indent amount for prettier print
    Returns
    ----------
    Metadata json
    """
    os.chmod(fl, 0o640)
    meta = load_json(fl)
    meta.update(data)
    with open(fl, 'wt') as fp:
        json.dump(meta, fp, indent=ind, sort_keys=True)
    os.chmod(fl, 0o440)

def main(subjs):
    for subj in subjs:
        for session in layout.get_sessions(): 
            fmaps = layout.get(subject=subj, datatype = 'fmap', session=session, extension='.json', return_type='file')
            if len(fmaps) > 0:
                process_fieldmaps(fmaps, subj, session)
            else:
                logging.warning('{} session {} does not have any fieldmap runs or this session does not exist for the subject.'.format(subj,session))     

# run across all bids subjects
# TODO: add configparser
print('starting')
subjs = layout.get_subjects()
main(subjs)
print('finished')
