"""
http://pep.databio.org/
"""

import os
import yaml
import itertools
import re

#import peppy
import pandas as pd


def config_info(config):
    single_cell = False
    subsamples = False
    multiple_projects = False
    multiple_flowcells = False
    
    for sample_id, sample_conf in config['samples'].items():
        if 'L00' in sample_conf.get('R1', ''):
            single_cell = True
        if ',' in sample_conf.get('R1', ''):
            subsamples = True
        if ',' in sample_conf.get('R2', ''):
            subsamples = True
        if ',' in sample_conf.get('Project_ID', ''):
            multiple_projects = True
        if ',' in sample_conf.get('Flowcell_ID', ''):
            multiple_flowcells = True
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
        drop_cols = []
        if info['multiple_flowcells']:
            drop_cols = ['Flowcell_Name', 'Flowcell_ID']
        if info['multiple_projects']:
            drop_cols.append('Project_ID')
        if drop_cols:
            df = df.drop(drop_cols, axis=1, errors='ignore')
    
    if 'R1' in df.columns:
        df['R1'] = 'R1'
    if 'R2' in df.columns:
        df['R2'] = 'R2'
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
        FC = sample_conf.get('Flowcell_Name', '').split(',')
        PID = sample_conf.get('Project_ID', '').split(',')
        run_number = 0
        for r1, r2, fc, pid in itertools.zip_longest(R1, R2, FC, PID):
            run_number += 1
            subsample_name = '{}_{}'.format(sample_id, run_number)
            subsamples[subsample_name] = {'subsample_name': subsample_name, 'sample_name': sample_id}
            if info['multiple_flowcells']:
                subsamples[subsample_name]['Flowcell_Name'] = fc
            if info['multiple_projects']:
                subsamples[subsample_name]['Project_ID'] = pid
            if info['single_cell']:
                patt = re.compile('(.*)\/(GCF-\d{4}-\d{3})\/(.*)\/.*_(S\d)_(L00\d)_R[1-2]_001.fastq.gz')
                m = patt.match(r1)
                if m:
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
    exp_params =  ['project_id', 'src_project_id', 'organism', 'machine', 'read_geometry']
    libprep_params = ['adapter', 'adapter2', 'read_orientation', 'workflow', 'libprep_name', 'delta_readlen']
    for k in exp_params:
        if k in config:
            exp_dict[k] = config[k]
    for n in libprep_params:
        if n in config:
            exp_dict[n] = config[n]
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
        sources['R2'] = '{Flowcell_Name}/{Project_ID}/{sample_name}/{sample_name}_{run_number}_{lane}_R1_001.fastq.gz'
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
        yaml.dump(peppy_conf, fh)    
        
