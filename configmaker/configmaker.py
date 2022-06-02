#!/usr/bin/env python

import sys
import os
import re
import glob
import logging
import json
import warnings
import subprocess
import re
import pprint
import shutil
from pathlib import Path
from collections.abc import Iterable
import six


import argparse
import pandas as pd
import oyaml as yaml

# enable local imports in script
#path_root = Path(__file__).parents[0]
#sys.path.append(str(path_root))

import descriptors 


SEQUENCERS = {
    'NB501038' : 'NextSeq 500',
    'SN7001334' : 'HiSeq 2500',
    'K00251' : 'HiSeq 4000',
    'M02675' : 'MiSeq NTNU',
    'M03942' : 'MiSeq StOlav',
    'M05617' : 'MiSeq SINTEF'
}

GCF_WORKFLOWS_SRC = "https://github.com/gcfntnu/gcf-workflows.git"

SNAKEFILE_TEMPLATE = """
from snakemake.utils import validate, min_version

pepfile:
    'pep/pep_config.yaml'
configfile:
    'config.yaml'

include:
    'src/gcf-workflows/{workflow}/{workflow}.smk'

"""

def setup_logger(verbose=False):
    logger = logging.getLogger('GCF-configmaker')
    fh = logging.FileHandler('.configmaker.debug')
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(levelname)s %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    if verbose:
        logger.setLevel(10)
        ch.setLevel(10)
        logger.debug('setting logging to debug ...')    
    return logger

logger = setup_logger()


def uniq_list():
    """
    list of uniq values order is preserved.
    """
    out = []
    [out.append(i) for i in seq if not out.count(i)]
    return out

class FullPaths(argparse.Action):
    """
    Expand user- and relative-paths.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        values = [os.path.abspath(os.path.expanduser(v)) for v in values]
        setattr(namespace, self.dest, values)
    
def is_dir(dirname):
    """
    Checks if a path is an actual directory.
    """
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_valid_gcf_id(arg, patt='GCF-\d{4}-\d{3}'):
    if arg is None:
        return True
    m = re.match(patt, arg)
    if m:
        return m.group().strip()
    else:
        msg = "{0} is not a valid GCF number (format: GCF-YYYY-NNN)".format(arg)
        raise argparse.ArgumentTypeError(msg)

def _match_project_dir(pth, project_id=None, test=False):
    """
    returns path and folder (project_id) of a project
    """
    if project_id:
        for fn in os.listdir(pth):
            if os.path.isdir(os.path.join(pth, fn)) and fn in project_id:
                return os.path.join(pth, fn), fn
        msg = "{0} is not present in run_folder: {1}".format(project_id, pth)
        logger.warning(msg)
        return None, None
    else:
        project_dir = None
        
        for fn in os.listdir(pth):
            
            if os.path.isdir(os.path.join(pth, fn)) and re.match('^GCF-\d{4}-\d{3}', fn):
                if project_dir is not None:
                    logger.error('runfolders contain more than one project folders existing: {}, other: {}'.format(project_id, fn) )
                    logger.error('Use `--project-id` option to choose one.')
                project_dir = os.path.join(pth, fn)
                project_id = fn
            elif test and re.match('GCF-\d{4}-\d{3}_samplesheet.tsv', fn):
                project_id = fn.split('_samplesheet.tsv')[0]
                project_dir = os.path.join(pth, project_id)
        if project_dir:
            logger.debug('project_dir match: {}'.format(project_dir))
            return project_dir, project_id
        raise RuntimeError('failed to identify any valid projects in runfolder: {}'.format(pth))



def get_data_from_samplesheet(fh):
    """
    returns [data] section as dataframe and [CustomOptions] section as key-value dict
    """
    custom_opts = False
    header = False
    header_d = {}
    opts_d = {}
    while True:
        line = fh.readline()
        if not line:
            msg = 'No [data]-section in samplesheet {}'.format(s.name)
            raise RuntimeError(msg)
        if line.startswith('[Data]'):
            return pd.read_csv(fh, dtype={'Sample_ID': str, 'Sample_Name': str}), opts_d, header_d
        elif line.startswith('[CustomOptions]'):
            custom_opts = True
            continue
        elif line.startswith('[Header]'):
            header = True
            continue
        elif custom_opts:
            els = [i.strip() for i in line.split(',') if i]
            if len(els) == 1:
                key, val = els[0], ''
            else:
                key, val = els[:2]
            if val.lower() == 'true':
                val = True
            if key == 'Organism' and val == 'N/A':
                val = None
            logger.debug('custom opt: {}:{}'.format(key, str(val)))
            opts_d[key] = val
        elif header:
            if line.startswith('['):
                header = False
                continue
            els = [i.strip() for i in line.split(',') if i]
            if len(els) == 1:
                key, val = els[0], ''
            else:
                key, val = els[:2]
            if val.lower() == 'true':
                val = True
            logger.debug('header: {}:{}'.format(key, str(val)))
            header_d[key] = val
            
def get_project_samples_from_samplesheet(args):
    """
    Return a dataframe containing project samples

    !! Assuming CustomOptions are equal between all samplesheets
    """
    
    samples_dataframe_list = []
    for samplesheet_i in args.samplesheet:
        with open(samplesheet_i, 'r') as fh:
            data, opts, header = get_data_from_samplesheet(fh)
            samples_dataframe_list.append(data)
    df = pd.concat(samples_dataframe_list)
    if args.project_id:
        # subset samples on project id
        keep = df.Sample_Project.isin(args.project_id)
        n_samples = df.shape[0]
        n_keep = sum(keep)
        if n_keep < n_samples:
            logger.debug('subsetting samples on project_id')
            logger.debug('{} samples kept out of {}'.format(n_keep, n_samples))
            df = df[df.Sample_Project.isin(args.project_id)]
    df['Sample_ID'] = df['Sample_ID'].astype(str)
    df['Project_ID'] = [[i] for i in df['Sample_Project']] # store project-ids as list
    df = df[['Sample_ID', 'Project_ID']]
    df = df.convert_dtypes()
    return df, opts, header

def match_fastq(sample_name, project_dir, rel_path=True):
    """
    Return fastq files matching a sample name.

    Returns paths relative to project directory
    """
    r1_fastq_files, r2_fastq_files = [],[]
    for fn in  os.listdir(project_dir):
        if fn == '{}_R1.fastq.gz'.format(sample_name):
            r1_fastq_files.extend([os.path.join(project_dir, fn)])
        elif fn == '{}_R2.fastq.gz'.format(sample_name):
            r2_fastq_files.extend([os.path.join(project_dir, fn)])
        elif fn == sample_name:
            r1_fastq_files.extend(glob.glob(os.path.join(project_dir, sample_name, sample_name + '*_R1_001.fastq.gz')))
            r2_fastq_files.extend(glob.glob(os.path.join(project_dir, sample_name, sample_name + '*_R2_001.fastq.gz')))
    if (len(r1_fastq_files) == 0) and (len(r2_fastq_files) == 0):
        warn_msg = 'Failed to match sample: {} with any fastq files in {}'.format(sample_name, project_dir)
        logger.warning(warn_msg)
        return None, None
    r1_fastq_files = sorted(r1_fastq_files)
    r2_fastq_files = sorted(r2_fastq_files)
    if rel_path:
        r1_fastq_files = [os.path.relpath(x, os.path.dirname(os.path.dirname(project_dir))) for x in r1_fastq_files]
        r2_fastq_files = [os.path.relpath(x, os.path.dirname(os.path.dirname(project_dir))) for x in r2_fastq_files]

    return r1_fastq_files, r2_fastq_files

def find_samples(df, args):
    """
    identify valid samples by existing fastq file names
    """
    sample_dict = {}
    project_dirs = [os.path.join(run_folder, project_id) for run_folder, project_id in zip(args.runfolders, args.project_id)]
    for index, row in df.iterrows():
        s_r1 = []
        s_r2 = []
        for p_pth in project_dirs:
            r1, r2 = match_fastq(row.Sample_ID, p_pth)
            if r1 is not None:
                s_r1.extend(r1)
            if r2 is not None:
                s_r2.extend(r2)
        if all([i is None for i in s_r1]) and all([i is None for i in s_r2]):
            warn_str = 'removing sample {} from SampleSheet due to missing fastq files!'.format(row.Sample_ID)
            logger.warning(warn_str)
        else:
            sample_dict[str(row.Sample_ID)] = {'R1': ','.join(s_r1),
                                               'R2': ','.join(s_r2),
                                               'Project_ID': ','.join(row.Project_ID),
                                               'Sample_ID': row.Sample_ID}
    return sample_dict


def find_samples_batch(df, project_dirs):
    """
    `find_samples` function adding Flowcell_ID postfix to Sample_ID
    """
    sample_dict = {}
    project_dirs = [os.path.join(run_folder, project_id) for run_folder, project_id in zip(args.runfolders, args.project_id)]
    for index, row in df.iterrows():
        for p_pth in project_dirs:
            r1, r2 = match_fastq(row.Sample_ID, p_pth)
            if (not r1) and (not r2):
                warn_str = 'sample {} not found in {}'.format(row.Sample_ID, p_pth)
                logger.warning(warn_str)
            else:
                r2 = [] if not r2 else r2
                Flowcell_ID = os.path.split(p_pth)[-2].split("_")[-1]
                Sample_ID = '{}_{}'.format(row.Sample_ID, Flowcell_ID)
                sample_dict[Sample_ID] = {'R1': ','.join(r1),
                                          'R2': ','.join(r2),
                                          'Project_ID': ','.join(row.Project_ID),
                                          'Sample_ID': Sample_ID,
                                          'Src_Sample_ID': row.Sample_ID}
    return sample_dict


def find_samples_test(df, project_dirs):
    sample_dict = {}
    for index, row in df.iterrows():
        sample_dict[str(row.Sample_ID)] = {'R1': '{}_R1.fastq.gz'.format(row.Sample_ID),
                                           'R2': '{}_R2.fastq.gz'.format(row.Sample_ID),
                                           'Project_ID': ','.join(row.Project_ID),
                                           'Sample_ID': row.Sample_ID}
    return sample_dict
                                           
def _customer_column_mapper(x):
    """
    map headers of customer-sheet to machine friendly headers
    """
    starts = [('Unique', 'Sample_ID'),
              ('External', 'External_ID'),
              ('Sample Group', 'Sample_Group'),
              ('Comment', 'Customer_Comments'),
              ('Sample biosource', 'Sample_Biosource'),
              ('Project', 'Project_ID'),
              ('Sample type', 'Sample_Type'),
              ('Sample Type', 'Sample_Type'),
              ('Index (If libraries are submitted  indicate what index sequence is used P7 )', 'Index1'),
              ('Index2', 'Index2'),
              ('Index1', 'Index1'),
              ('Sequence1', 'Index_Sequence1'),
              ('Sequence2', 'Index_Sequence2'),
              ('Plate location', 'Plate'),
              ('Sample buffer', 'Sample_Buffer'),
              ('Sample Buffer', 'Sample_Buffer'),
              ('Volume', 'Volume'),
              ('Quantification', 'Quantification'),
              ('Concentration', 'Concentration'),
              ('260/280', '260/280'),
              ('260/230', '260/230'),
              ('Organism', 'Organism'),
              ('RIN', 'RIN')
              ]
    for src, dst in starts:
        if x.startswith(src):
            return dst
    # unknown header value (may be customer added)
    src_sanitized = x.title()
    remove = """- ? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .""".split()
    for r in remove:
        src_sanitized = src_sanitized.replace(r, '')
        src_sanitized = src_sanitized.replace(' ', '_')
    return 'Submitted_' + src_sanitized

def _lab_column_mapper(x):
    """
    map headers of wetlab-sheet to machine friendly headers
    """
    starts = [('Concentration', 'Concentration'),
              ('260/280', '260/280'),
              ('260/230', '260/230'),
              ('Comment', 'Comments'),
              ('Sample_ID', 'Sample_ID'),
              ('Project', 'Project_ID'),
              ('RIN', 'RIN'),
              ('SpikeIn', 'SpikeIn'),
              ('Fragment_Length', 'Fragment_Length'),
              ('Fragment_SD', 'Fragment_SD'),
              ('Sample_Name', 'Lab_Sample_Name'),
              ('KIT', 'KIT'),
              ('ERCC', 'ERCC')
              ]
    for src, dst in starts:
        if x.startswith(src):
            return dst
    src_sanitized = x.title()
    remove = """- ? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .""".split()
    for r in remove:
        src_sanitized = src_sanitized.replace(r, '')
        src_sanitized = src_sanitized.replace(' ', '_')
    return 'Lab_' + src_sanitized

def read_customer_sheet(fn):
    _dtypes = {'Unique Sample ID': str, 'External ID (optional reference sample ID)': str}
    df = pd.read_excel(fn, sheet_name='Sample-Submission-Form', skiprows=14, dtype=_dtypes)
    desc = descriptors.descriptors.findall_header_descriptors(df, mapper=_customer_column_mapper) # identify any header descriptors
    df = df.rename(columns=_customer_column_mapper)
    remove_cols = ['Sample_Type', 'Sample_Buffer', 'Volume', 'Quantification'] # wetlab only
    remove_cols = list(set(df.columns).intersection(remove_cols))
    df = df.drop(remove_cols, axis=1)
    df = df.replace('NA', pd.NA)
    df = df.dropna(axis='columns', how='all') # remove empty cols
    if not df.empty:
        df = df.convert_dtypes()
    logger.debug('customer descriptors: {}'.format(str(desc)))
    return df, desc

def read_lab_sheet(fn):
    df = pd.read_excel(fn, sheet_name='INFO (GCF-lab only)', dtype={'Sample_ID': str})
    desc = descriptors.descriptors.findall_header_descriptors(df, mapper=_lab_column_mapper)
    df = df.rename(columns=_lab_column_mapper)
    legacy_cols = list(set(['Sample_Name','KIT']).intersection(df.columns))
    df = df.drop(legacy_cols, axis=1, errors='ignore')
    df = df.replace('NA', pd.NA)
    df = df.dropna(axis='columns', how='all') # remove empty cols
    if not df.empty:
        df = df.convert_dtypes()
    logger.debug('lab descriptors: {}'.format(str(desc)))
    return df, desc


def sample_submission_form_parser(ssub_path, keep_batch=None):
    """read submission form excel file

    merge sheets (customer + lab) and santize column names, values and add column descriptors

    returns a merged dataframe from lab and customer where each row represent a sample and a descriptor dictionary keyed in column names
    """

    
    customer, desc = read_customer_sheet(ssub_path)
    lab, lab_desc = read_lab_sheet(ssub_path)
    # lab-sheet will take presedence over customer filled columns
    shared_cols = list(set(customer.columns).intersection(lab.columns))
    if 'Sample_ID' in shared_cols:
        shared_cols.remove('Sample_ID')
    customer = customer.drop(shared_cols, axis=1)
    for k, v in lab_desc.items():
        if v:
            if k not in shared_cols and not desc.get(k):
                logger.debug('lab descriptor: {}:{}'.format(k, v))
                desc[k] = v
    flowcell_name = os.path.basename(os.path.split(ssub_path)[0])
    flowcell_id = flowcell_name.split('_')[-1]
    customer['Flowcell_Name'] = flowcell_name #flowcell folder name
    customer['Flowcell_ID'] = flowcell_id

    if keep_batch:
        customer['Sample_ID'] = customer['Sample_ID'].astype(str) + "_" + flowcell_id
        lab['Sample_ID'] = lab['Sample_ID'].astype(str) + "_" + flowcell_id
        
    if not lab.empty:
        merged_ssub = pd.merge(customer, lab, on='Sample_ID', how='inner')
    else:
        merged_ssub = customer
    merged_ssub["Sample_ID"] = merged_ssub["Sample_ID"].astype(str)
    merged_ssub.index = merged_ssub["Sample_ID"]
    desc = descriptors.descriptors.add_default_descriptors(merged_ssub, desc)
 
    return merged_ssub, desc

def _make_header_uniq(df):
    columns = df.columns.values.copy()
    dup_vals = set(df.columns[df.columns.duplicated()])
    logger.error('duplicate col names: {}'.format(dup_vals))
    for d in dup_vals:
        dups = df.columns[df.columns==d]
        uniq_names = ['{}_{}'.format(k, i+1) for i, k in enumerate(dups)]
        columns[df.columns==d] = uniq_names
    df.columns = columns
    return df

def merge_samples_with_submission_form(sample_dict, args):
    """
    """
    submission_forms, desc = {}, {}
    for submission_form in args.ssub:
        ssub, sub_desc = sample_submission_form_parser(submission_form, keep_batch=args.keep_batch)
        pth = os.path.abspath(submission_form)
        if len(set(ssub.columns)) != len(ssub.columns):
            ssub = _make_header_uniq(ssub)
        submission_forms[pth] = ssub.to_dict(orient='index')

        for k, v in sub_desc.items():
            if k in desc:
                if desc[k] != v:
                    default_v = descriptors.descriptors.DEFAULT_DESCRIPTORS.get(k)
                    if desc[k] == default_v:
                        # update descriptor if is new and not default value
                        desc[k] = v
            else:
                desc[k] = v

    merge = dict()
    for pth, sf_dict in submission_forms.items():
        for sample_id, vals in sf_dict.items():
            if sample_id not in merge:
                merge[sample_id] = vals
            else:
                # merge info from sample
                sample = merge[sample_id]
                for k, v in vals.items():
                    if k in ['Flowcell_Name', 'Flowcell_ID', 'Project_ID']:
                        # special case flowcell name to multiple values by comma sep 
                        v = ','.join([sample[k], v])
                    elif sample[k] == v or (pd.isnull(v) and pd.isnull(sample[k])):
                        # equal info between submission forms
                        pass
                    else:
                        # unequal info in submission forms and we are not looking at flowcell
                        msg = "Sampleinfo ({}) at {} are updated with values from {}/Sample-Submission-Form.xlsx. ".format(v, k, pth)
                        msg2 = "Specify a custom sample submission form with --sample-submission-form to force values."
                        logger.warning(msg + msg2)
                    sample[k] = v
                merge[sample_id] = sample
    merge = pd.DataFrame.from_dict(merge, orient='index')
    check_existence_of_samples(sample_dict.keys(), merge)
    sample_df = pd.DataFrame.from_dict(sample_dict, orient='index')
    if 'Project_ID' in sample_df.columns and 'Project_ID' in merge:
        # use Project_ID from samplesheet over sample-submission-form
        merge = merge.drop('Project_ID', axis=1)
    sample_df = sample_df.merge(merge, on='Sample_ID', how='left')
    sample_df.reset_index()
    sample_df.index = sample_df['Sample_ID']
        
    if args.new_project_id:
        sample_df = sample_df.rename(columns={'Project_ID': 'Src_Project_ID'})
        sample_df['Project_ID'] = args.new_project_id
    
    desc = descriptors.descriptors.add_default_descriptors(sample_df, desc)
    sample_df, desc = descriptors.descriptors.infer_by_descriptor(sample_df, desc)
    if 'Organism' in sample_df.columns:
        if (len(set(sample_df.Organism.values)) == 1): # single customer org
            if args.organism is not None:
                sample_df['Organism'] = args.organism
            else:
                # if org is N/A in samplesheet but single org in sample_df
                args.organism = list(sample_df.Organism)[0]
    return sample_df, desc

def check_existence_of_samples(samples, df):
    diff = list(set(samples) - set(df['Sample_ID'].astype(str)))
    if diff:
        extra = ''
        n_diff = len(diff)
        if n_diff > 10:
            diff = diff[:10]
            extra = ',... +{} samples'.format(n_diff - 10)
        vals = ','.join(list(diff)) + extra
        
        logger.warning("Samples {} are contained in SampleSheet, but not in sample submission form. Sample info from these samples will have empty values.".format(vals))
    diff2 = list(set(df['Sample_ID'].astype(str)) - set(samples))                   
    if diff2:
        extra = ''
        n_diff2 = len(diff2)
        if n_diff2 > 10:
            diff2 = diff2[:10]
            extra = ',... +{} samples'.format(n_diff2 - 10)
        vals2 = ','.join(list(diff2)) + extra
        logger.warning("Samples {} are contained in sample submission form, but not in SampleSheet. Sample info from these samples are omitted.".format(vals2))
    return None

def find_read_geometry(runfolders):
    n_matches = set()
    for fn in runfolders:
        stats_fn = os.path.join(fn, 'Stats', 'Stats.json')
        read_geometry = []
        with open(stats_fn) as fh:
            S = json.load(fh)
        for read in S['ReadInfosForLanes'][0]['ReadInfos']:
            if not read['IsIndexedRead']:
                n_cycles = int(read['NumCycles'])
                read_geometry.append(n_cycles)
        n_matches.add(','.join(map(str, read_geometry)))
    if len(n_matches) > 1:
        raise ValueError('Read geometry mismatch between runfolders. Check Stats.json!')
    return read_geometry

def find_machine(runfolders):
    n_matches = set()
    for pth in runfolders:
        machine_code = os.path.basename(pth).split('_')[1]
        machine = SEQUENCERS.get(machine_code, '')
        n_matches.add(machine)
    if len(n_matches) > 1:
        #logger.warning('Multiple sequencing machines identified!')
        raise ValueError('Multiple sequencing machines identified!')
    return machine

def create_default_config(merged_samples, opts, args, fastq_dir=None, descriptors=None, write_yaml=True):
    """
    create configuration dictionary
    """
    samples = merged_samples.to_dict(orient='index')
    config = {}
    if len(set(args.project_id)) == 1:
        args.project_id = [args.project_id[0]]
    
    if args.new_project_id:
         config['project_id'] = args.new_project_id
         config['src_project_id'] = args.project_id
    else:
         config['project_id'] = args.project_id

    if 'Organism' in opts:
        config['organism'] = opts['Organism']
    if args.organism is not None:
        if pd.isnull(args.organism) or args.organism in ['N/A', 'NA', '<NA>', '', None]:
            config['organism'] = 'N/A'
        else:
            config['organism'] = args.organism

    if 'Libprep' in opts:
        config['libprepkit'] = opts['Libprep']
    if args.libkit is not None:
        config['libprepkit'] = args.libkit

    config['read_geometry'] = find_read_geometry(args.runfolders)
    config['machine'] = args.machine or find_machine(args.runfolders)

    batch = {}
    if args.keep_batch:
        batch['name'] = 'Flowcell_ID'
    else:
        batch['method'] = 'skip'
    quant = {'batch': batch}
    config['quant'] = quant

    if fastq_dir:
        config['fastq_dir'] = fastq_dir

    if descriptors:
        current_col_names = set()
        for k, v in samples.items():
            for n in v.keys():
                current_col_names.add(n)
        keys = current_col_names.intersection(descriptors.keys())
        config['descriptors'] = {}
        for k in keys:
            if k in descriptors:
                if descriptors[k]:
                    config['descriptors'][k] = descriptors[k]
        
    config['samples'] = {}
    for sample_id, col in samples.items():
        config['samples'][sample_id] = {}
        for col_name, val in col.items():
            if isinstance(val, Iterable) and not isinstance(val, six.string_types):
                val = list(map(str, val))
            else:
                val = str(val)
            config['samples'][sample_id][col_name] = val

    if write_yaml:
        yaml.dump(config, args.output, default_flow_style=False, tags=False)
    
    return config

def create_fastq_dir(sample_dict, args, output_dir=None, overwrite=True):
    """
    Symlink fastq files into project specific names and destinations
    """
    if not args.create_fastq_dir:
        return None
    default_fastq_dir = output_dir or os.path.join('data', 'raw', 'fastq')
    if os.path.exists(default_fastq_dir) and overwrite:
        shutil.rmtree(default_fastq_dir)
    os.makedirs(default_fastq_dir, exist_ok=True)

    s_ids = sample_dict.keys()
    if args.keep_batch:
        # split out date addition postfix from sample_ids with support for sample_ids with underscore in name
        s_ids = ['_'.join(n[:-1]) if len(n)>2 else n[0] for n in [e.split('_') for e in s_ids]]
        s_ids = list(set(s_ids)) # unique sample_ids

    project_dirs = [os.path.join(run_folder, project_id) for run_folder, project_id in zip(args.runfolders, args.project_id)]
    for sample_id in s_ids:
        for pid in project_dirs:
            r1_src, r2_src = match_fastq(sample_id, pid, rel_path=False)
            r1_dst, r2_dst = match_fastq(sample_id, pid, rel_path=True)
            if not any([r1_src, r1_dst, r2_src, r2_dst]):
                continue
            for src, dst in zip(r1_src, r1_dst):
                if src is not None:
                    dst = os.path.join(default_fastq_dir, dst)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    os.symlink(src, dst)
            for src, dst in zip(r2_src, r2_dst):
                if src is not None:
                    dst = os.path.join(default_fastq_dir, dst)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    os.symlink(src, dst)
    return default_fastq_dir    

def create_project(config, src_dir=None):
    """
    create GCF project structure
    """
    src_dir = src_dir or 'src'
    wf_path = os.path.join(src_dir, 'gcf-workflows')
    if not os.path.exists(src_dir):
        os.makedirs(src_dir, exist_ok=True)
        cmd = 'cd src && git clone {}'.format(GCF_WORKFLOWS_SRC)
        subprocess.check_call(cmd, shell=True)

    with open(os.path.join(wf_path, 'libprep.config'), 'r') as libprepconf_fh:
        libconf = yaml.load(libprepconf_fh, Loader=yaml.FullLoader)

    libkit = config['libprepkit'] + (' PE' if len(config['read_geometry']) > 1 else ' SE')
    kitconf = libconf.get(libkit, None)
    if not kitconf:
        logger.warning('Libprepkit {} is not defined in libprep.config. Running with default settings.'.format(libkit))
        workflow = 'default'
    else:
        workflow = kitconf['workflow']

    if not 'workflow' in config:
        config['workflow'] = workflow
        
    with open('Snakefile', 'w') as sn:
        sn.write(SNAKEFILE_TEMPLATE.format(workflow=workflow))

def project_summary(config):
    """
    print project summary
    """
    summary = dict()
    for s, info in config['samples'].items():
        summary[s] = set(x.split("/")[0] for x in info['R1'].split(","))

    count = dict()
    for s,f in summary.items():
        if len(f) in count.keys():
            count[len(f)] += 1
        else:
            count[len(f)] = 1
    dirname = os.path.dirname(str(args.output))
    with open(os.path.join(dirname, ".configmaker.log"), "w") as conflog:
        print("Summary:")
        for nf, ns in count.items():
            line = "{} sample{} found in {} flowcell{}".format(ns, "s" if ns > 1 else "", nf, "s" if nf > 1 else "")
            print(line)
            conflog.write(line + "\n")
        conflog.write("Sample summary:\n")
        for s, f in summary.items():
            conflog.write("Sample {} found in: {}\n".format(s, ', '.join(f)))
    print("Full sample summary log written to {}".format(os.path.join(dirname, ".configmaker.log")))



def check_input(args):
    logger.debug('running check_input ...')
    dirs, ids, samplesheets, submission_forms = [],[],[],[]
    for pth in args.runfolders:
        project_dir, project_id = _match_project_dir(pth, project_id=args.project_id, test=args.test)
        if project_id:
            dirs.append(project_dir)
            ids.append(project_id)
            samplesheet_fn = os.path.join(pth, 'SampleSheet.csv')
            ssub_fn = os.path.join(pth, 'Sample-Submission-Form.xlsx')
            if os.path.exists(samplesheet_fn) and args.samplesheet is None:
                logger.debug('identified samplesheet: {}'.format(samplesheet_fn))
                samplesheets.append(samplesheet_fn)
            if os.path.exists(ssub_fn) and args.ssub is None:
                logger.debug('identified submission form: {}'.format(ssub_fn))
                submission_forms.append(ssub_fn)
            
    if args.samplesheet is not None:
        logger.debug('overriding samplesheet with command line arg: {}'.format(args.samplesheet.name))
        args.samplesheet = [args.samplesheet.name]
    else:
        if len(samplesheets) == 0:
            msg = "cannot find SampleSheet.csv in runfolders. Use --samplesheet for manual override"
            logger.error(msg)
            raise RuntimeError(msg) 
        args.samplesheet = samplesheets
    if args.ssub is not None:
        args.ssub = [args.ssub.name]
    else:
        if len(submission_forms) == 0:
            msg = "cannot find Sample-Submission-Form.xlsx in runfolders. Use --submission-form for manual override"
            logger.error(msg)
            raise RuntimeError(msg) 
        args.ssub = submission_forms
    args.project_id = ids
    args.runfolders = [os.path.dirname(p) for p in dirs]

    args.organism = descriptors.fuzzmatch.fuzzmatch_organism(args.organism)
    
    return args
            


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--project-id", nargs="+", help="Project ID", default=None, type=is_valid_gcf_id)
    parser.add_argument("-P", "--new-project-id", help="New Project ID", default=None, type=is_valid_gcf_id)
    parser.add_argument("runfolders", nargs="+", help="Path(s) to flowcell dir(s)", action=FullPaths, type=is_dir)
    parser.add_argument("-s", "--sample-sheet", dest="samplesheet", type=argparse.FileType('r'), help="IEM Samplesheet")
    parser.add_argument("-o", "--output", default="config.yaml", help="Output config file", type=argparse.FileType('w'))
    parser.add_argument("-S", "--sample-submission-form", dest="ssub", type=argparse.FileType('r'), help="GCF Sample Submission Form")
    parser.add_argument("--organism",  help="Organism (if applicable to all samples). Overrides value from samplesheet.")
    parser.add_argument("--libkit",  help="Library preparation kit name. (if applicable for all samples). Overrides value from samplesheet.")
    parser.add_argument("--machine",  help="Sequencer model.")
    parser.add_argument("--create-fastq-dir", action='store_true', help="Create fastq dir and symlink fastq files")
    parser.add_argument("--create-project", action='store_true', help="Pull analysis pipeline and snakemake file based on libkit")
    parser.add_argument("--skip-peppy", action='store_true', help="Skip creation of a peppy project")
    parser.add_argument("--keep-batch", action='store_true', help="Sample names will be made unique for each batch.")
    parser.add_argument("--verbose", action='store_true', help="Verbose output")
    parser.add_argument("--test", action='store_true', help="Activate test-mode. (no fastq files needed)")

    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = parse_args()
    if args.verbose:
        logger = setup_logger(verbose=True)
    args = check_input(args)
    samples_df, custom_opts, header = get_project_samples_from_samplesheet(args)
    args.organism = args.organism or custom_opts.get('Organism')
    args.organism = args.organism or descriptors.fuzzmatch.fuzzmatch_organism(args.organism)
        

    if args.test:
        sample_dict = find_samples_test(samples_df, args)
    else:
        if args.keep_batch:
            sample_dict = find_samples_batch(samples_df, args)
        else:
            sample_dict = find_samples(samples_df, args)

    merged_samples, desc = merge_samples_with_submission_form(sample_dict, args)

    _ = pd.DataFrame.from_dict(desc).T
    summary = merged_samples.dtypes.rename('pd_dtype').to_frame().merge(_, how='left', left_index=True, right_index=True)
    print(summary.dropna(axis='columns', how='all').fillna(''))


    fastq_dir = create_fastq_dir(sample_dict, args)

    config = create_default_config(merged_samples, custom_opts, args, fastq_dir=fastq_dir, descriptors=desc)

    if args.create_project:
        create_project(config)

    if not args.skip_peppy:
        import peppy_support
        peppy_support.create_peppy(config, output_dir='pep')
    

    project_summary(config)
