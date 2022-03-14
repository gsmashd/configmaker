"""
http://pep.databio.org/
"""

import os
import yaml

#import peppy
import pandas as pd

def _contains_subsamples(config):
    """returns true if project contains subsamples

    mutiple fastq files of same readtype within one sample are subsamples
    """
    for sample_id, sample_conf in config['samples'].items():
        if ',' in sample_conf.get('R1', ''):
            return True
        if ',' in sample_conf.get('R2', ''):
            return True
    return False


def conifg2sampletable(config, subsamples=False, drop_empty_cols=True):
    """returns sampletable dataframe
    """
    samples = config['samples']
    df = pd.DataFrame.from_dict(samples, orient='index')
    if drop_empty_cols:
        # drop columns with all NaN or empty string
        df = df.dropna(axis='columns', how='all')
        df = df.loc[:,(df == '').sum(0) < df.shape[0]]
        
    if subsamples:
        df = df.drop(columns={'Flowcell_Name'})
    
    if 'R1' in df.columns:
        df['R1'] = 'R1'
    if 'R2' in df.columns:
        df['R2'] = 'R2'
    df = df.set_index('Sample_ID')
    df = df.reset_index()
    df = df.rename(columns={'Sample_ID': 'sample_name'})
    
    return df

def config2subsampletable(config):
    """returns subsampletable dataframe
    """
    subsamples = {}
    for sample_id, sample_conf in config['samples'].items():
        subsample = {}
        R1 = sample_conf.get('R1', '').split(',')
        R2 = sample_conf.get('R2', '').split(',')
        FC = sample_conf.get('Flowcell_Name', '').split(',')
        run_number = 0
        for r1, r2, fc in zip(R1, R2, FC):
            run_number += 1
            subsample_name = '{}_{}'.format(sample_id, run_number)
            subsamples[subsample_name] = {'Flowcell_Name': fc, 'subsample_name': subsample_name, 'sample_name': sample_id}
    if len(subsamples) > 0:
        return pd.DataFrame.from_dict(subsamples, orient='index')[['sample_name', 'subsample_name', 'Flowcell_Name']]
    else:
        return None

def peppy_project_dict(config, contains_subsamples=False):
    """return project config dictionary
    """
    
    peppy_project = {}
    peppy_project['pep_version'] = '2.0.0'
    peppy_project['sample_table'] = 'sample_table.csv'
    if contains_subsamples:
        peppy_project['subsample_table'] = 'subsample_table.csv'
    derive = {}
    derive['attributes'] = ['R1', 'R2']
    sources = {}
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
    subsamples = _contains_subsamples(config)
    sampletable = conifg2sampletable(config, subsamples)
    subsampletable = config2subsampletable(config)

    peppy_conf = peppy_project_dict(config, contains_subsamples=subsamples)

    # write output
    os.makedirs(peppy_dir, exist_ok=True)
    sampletable.to_csv(os.path.join(peppy_dir, 'sample_table.csv'), index=None)
    if subsamples:
        subsampletable.to_csv(os.path.join(peppy_dir, 'subsample_table.csv'), index=None)
    with open(os.path.join(peppy_dir, 'project_config.yaml'), 'w') as fh:
        yaml.dump(peppy_conf, fh)    
        
