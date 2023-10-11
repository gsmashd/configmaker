"""
http://pep.databio.org/
"""

import os
import oyaml as yaml
import itertools
import re
import sys

#import peppy
import pandas as pd

BASE = ['project_id',
        'src_project_id',
        'organism']

RUN =  ['workflow',
        'machine',
        'read_geometry']

LIBPREP = ['adapter',
           'adapter2',
           'read_orientation',
           'libprepkit',
           'molecule',
           'library_strategy',
           'library_selection',
           'library_source',
           'library_strand']

STUDY = ['experiment_title',
         'experiment_summary',
         'experiment_principal_inverstigator',
         'experiment_contributor']



def config_info(config):
    single_cell = False
    subsamples = False
    multiple_projects = False
    multiple_flowcells = False
    
    for sample_id, sample_conf in config['samples'].items():
        if 'I1' in sample_conf:
            single_cell = True
        if ',' in sample_conf.get('R1', ''):
            subsamples = True
        if ',' in sample_conf.get('R2', ''):
            subsamples = True
        if ',' in sample_conf.get('Project_ID', ''):
            multiple_projects = True
        if ',' in sample_conf.get('Flowcell_ID', ''):
            multiple_flowcells = True
            #if len(sample_conf.get('Flowcell_ID').split(',')) == len(sample_conf.get('R1').split(',')):
            #    subsamples = False
    return {'subsamples': subsamples, 'multiple_projects':multiple_projects, 'multiple_flowcells': multiple_flowcells, 'single_cell':single_cell}

def _empty_col(col):
    if all(col==''):
        return True
    if col.isna().all():
        return True
    if col.isnull().all():
        return True
    return False

def conifg2sampletable(config, info, drop_empty_cols=True):
    """returns sampletable dataframe
    """
    samples = config['samples']
    df = pd.DataFrame.from_dict(samples, orient='index')
    if drop_empty_cols:
        # drop columns with all NaN or empty string
        empty = df.apply(_empty_col, axis=0) # rm empty cols
        df = df.loc[:,~empty]
        
    if info['subsamples']:
        drop_cols = [i for i in df.columns if i.endswith('_md5sum')]
        #drop_cols = []
        if info['multiple_flowcells']:
            drop_cols.extend(['Flowcell_Name', 'Flowcell_ID'])
        if info['multiple_projects']:
            drop_cols.append('Project_ID')
        if drop_cols:
            df = df.drop(drop_cols, axis=1, errors='ignore')
    
    if 'R1' in df.columns:
        df['R1'] = 'R1'
    if 'R2' in df.columns:
        df['R2'] = 'R2'
    if 'I1' in df.columns:
        df['I1'] = 'I1'
    df = df.set_index('Sample_ID')
    df = df.reset_index()
    df = df.rename(columns={'Sample_ID': 'sample_name'})
    
    return df

def config2subsampletable(config, info):
    """returns subsampletable dataframe
    """
    subsamples = {}
    for sample_id, sample_conf in config['samples'].items():
        subsample = {}
        R1 = sample_conf.get('R1', '').split(',')
        R2 = sample_conf.get('R2', '').split(',')
        I1 = sample_conf.get('I1', '').split(',')
        R1_MD5 = sample_conf.get('R1_md5sum', '').split(',')
        R2_MD5 = sample_conf.get('R2_md5sum', '').split(',')
        I1_MD5 = sample_conf.get('I1_md5sum', '').split(',')
        FC = sample_conf.get('Flowcell_Name', '').split(',')
        PID = sample_conf.get('Project_ID', '').split(',')
        run_number = 0
        for r1, fc, pid in itertools.zip_longest(R1, FC, PID):
            run_number += 1
            subsample_name = '{}_{}'.format(sample_id, run_number)
            subsamples[subsample_name] = {'subsample_name': subsample_name, 'sample_name': sample_id}
            if info['subsamples']:
                if len(R1_MD5) == len(R1): 
                    subsamples[subsample_name]['R1_md5sum'] = R1_MD5[run_number-1]
                if len(R2_MD5) == len(R2) and R2[0]:
                    subsamples[subsample_name]['R2_md5sum'] = R2_MD5[run_number-1]
                if len(I1_MD5) == len(I1) and I1[0]:
                    subsamples[subsample_name]['I1_md5sum'] = I1_MD5[run_number-1]
            
            if info['multiple_flowcells']:
                subsamples[subsample_name]['Flowcell_Name'] = fc
            if info['multiple_projects']:
                subsamples[subsample_name]['Project_ID'] = pid
            if info['single_cell']:
                m = re.match('(.*)\/(GCF-\d{4}-\d{3})\/(.*)\/.*_(S\d+)_(L00\d)_R[1-2]_001.fastq.gz', r1)
                if m:
                    _fc, _pid, _sample_id, run_id, lane = m.groups()
                else:
                    m = re.match('(.*)\/(GCF-\d{4}-\d{3})\/(.*)_(S\d+)_(L00\d)_R[1-2]_001.fastq.gz', r1)
                    _fc, _pid, _sample_id, run_id, lane = m.groups()
                subsamples[subsample_name]['lane'] = lane
                subsamples[subsample_name]['run_number'] = run_id
                
            
    if len(subsamples) > 0:
        df = pd.DataFrame.from_dict(subsamples, orient='index')
        df = df.set_index('sample_name').reset_index()
        return df
    else:
        return None

def config2experimentinfo(config):
    """extract global experiment/data info
    """
    exp_dict = {}
    for k in BASE:
        if k in config:
            exp_dict[k] = config[k]
    for n in STUDY:
        if n in config:
            exp_dict[n] = config[n]        
    for k in RUN:
        if k in config:
            exp_dict[k] = config[k]
    for n in LIBPREP:
        if n in config:
            exp_dict[n] = config[n]
    
    exp_dict['protocol'] = config.get('protocol', {})
    exp_dict['descriptors'] = config.get('descriptors', {})
    
    return exp_dict
        

def peppy_project_dict(config, info):
    """return project config dictionary
    """
    
    peppy_project = {}
    peppy_project['pep_version'] = '2.0.0'
    peppy_project['sample_table'] = 'sample_table.csv'
    global_params = config2experimentinfo(config)
    peppy_project.update(global_params)
    
    if info['subsamples']:
        peppy_project['subsample_table'] = 'subsample_table.csv'
    derive = {}
    derive['attributes'] = ['R1', 'R2']
    sources = {}
    if info['single_cell']:
        # bcl2fastq without --no-lane-splitting
        sources['R1'] = '{Flowcell_Name}/{Project_ID}/{sample_name}/{sample_name}_{run_number}_{lane}_R1_001.fastq.gz'
        sources['R2'] = '{Flowcell_Name}/{Project_ID}/{sample_name}/{sample_name}_{run_number}_{lane}_R2_001.fastq.gz'
        sources['I1'] = '{Flowcell_Name}/{Project_ID}/{sample_name}/{sample_name}_{run_number}_{lane}_I1_001.fastq.gz'
    else:
        sources['R1'] = '{Flowcell_Name}/{Project_ID}/{sample_name}_R1.fastq.gz'
        sources['R2'] = '{Flowcell_Name}/{Project_ID}/{sample_name}_R2.fastq.gz'
    derive['sources'] = sources
    peppy_project['sample_modifiers'] = {}
    peppy_project['sample_modifiers']['derive'] = derive
    return peppy_project
    

def create_peppy(config, output_dir='peppy_project'):
    """createa a peppy project directory
    """
    peppy_project = {}
    peppy_dir = output_dir or 'peppy_project'
    info = config_info(config)
    sampletable = conifg2sampletable(config, info)
    subsampletable = config2subsampletable(config, info)
    peppy_conf = peppy_project_dict(config, info)

    # write output
    os.makedirs(peppy_dir, exist_ok=True)
    sampletable.to_csv(os.path.join(peppy_dir, 'sample_table.csv'), index=None)
    if info['subsamples']:
        subsampletable.to_csv(os.path.join(peppy_dir, 'subsample_table.csv'), index=None)
    with open(os.path.join(peppy_dir, 'pep_config.yaml'), 'w') as fh:
        yaml.safe_dump(peppy_conf, fh)    
        
