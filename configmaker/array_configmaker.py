#!/usr/bin/env python

import sys
import os
import re
import glob
import argparse
import pandas as pd
import yaml
import logging
import json
import warnings
import subprocess

logger = logging.getLogger('GCF-configmaker')
logger.setLevel(logging.WARNING)

GCF_WORKFLOWS_SRC = "https://github.com/gcfntnu/gcf-workflows.git"

SNAKEFILE_TEMPLATE = """
from snakemake.utils import validate, min_version

configfile:
    'config.yaml'

include:
    'src/gcf-workflows/{workflow}/{workflow}.smk'

"""

from . import configmaker

            
def get_project_samples_from_samplesheet(samplesheet, runfolders, project_id):
    """
    Return a dataframe containing project samples
    """
    ss = configmaker.inspect_samplesheet(samplesheet, runfolders)
    df_list = []
    for sheet in ss:
        with open(sheet, 'r') as s:
            data, opts = configmaker.get_data_from_samplesheet(s)
            df_list.append(data)
    df = pd.concat(df_list)
    if 'Sample_Group' in df.columns:
        df.drop('Sample_Group', axis='columns', inplace=True)
    if project_id:
        if 'Sample_Project' in df.columns:
            df = df[df.Sample_Project.isin(project_id)]
        else:
            df['Sample_Project'] = project_id or ''
    
    df['Sample_ID'] = df['Sample_ID'].astype(str)
    df = df[['Sample_ID']]
    df = df.drop_duplicates(['Sample_ID'])
    return df, opts

def create_default_config(sample_dict, opts, args, idat_dir=None):
    config = {}

    if args.new_project_id:
         config['project_id'] = args.new_project_id
         config['src_project_id'] = args.project_id
    else:
         config['project_id'] = args.project_id

    if 'Organism' in opts:
        config['organism'] = opts['Organism']
    if args.organism is not None:
        config['organism'] = args.organism

    if 'Libprep' in opts:
        config['libprepkit'] = opts['Libprep']
    if args.libkit is not None:
        config['libprepkit'] = args.libkit

    config['machine'] = args.machine

    batch = {}
    if args.keep_batch:
        batch['name'] = 'Plate'
    else:
        batch['method'] = 'skip'
    quant = {'batch': batch}
    config['quant'] = quant

    if idat_dir:
        config['fastq_dir'] = config['idat_dir'] = idat_dir

    config['samples'] = sample_dict

    return config

def match_idat():
    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--project-id", nargs="+", help="Project ID", default=None, type=is_valid_gcf_id)
    parser.add_argument("-P", "--new-project-id", help="New Project ID", default=None, type=is_valid_gcf_id)
    parser.add_argument("runfolders", nargs="+", help="Path(s) to flowcell dir(s)", action=FullPaths, type=is_dir)
    parser.add_argument("-s", "--sample-sheet", dest="samplesheet", type=argparse.FileType('r'), help="IEM Samplesheet")
    parser.add_argument("-o", "--output", default="config.yaml", help="Output config file", type=argparse.FileType('w'))
    parser.add_argument("-S", "--sample-submission-form", dest="ssub", type=argparse.FileType('r'), help="GCF Sample Submission Form")
    parser.add_argument("--organism",  help="Organism (if applicable to all samples). Overrides value from samplesheet.", default='homo_sapiens')
    parser.add_argument("--libkit",  help="Library preparation kit name. (if applicable for all samples). Overrides value from samplesheet.", defualt='epic')
    parser.add_argument("--machine",  help="Sequencer model.", default="iscan")
    parser.add_argument("--create-idat-dir", action='store_true', help="Create idat dir and symlink files")
    parser.add_argument("--create-project", action='store_true', help="Pull analysis pipeline and snakemake file based on libkit")
    parser.add_argument("--keep-batch", action='store_true', help="Sample names will be made unique for each batch.")

    args = parser.parse_args()

    project_dirs, args.project_id = configmaker.inspect_dirs(args.runfolders, args.project_id)
    s_df, opts = get_project_samples_from_samplesheet(args.samplesheet, args.runfolders, args.project_id)
    if args.keep_batch:
        sample_dict = configmaker.find_samples_batch(s_df, project_dirs)
    else:
        sample_dict = configmaker.find_samples(s_df, project_dirs)

    ssub_d = dict()
    if args.ssub is None:
        for pth in args.runfolders:
            ssub_fn = os.path.join(pth, 'Sample-Submission-Form.xlsx')
            if os.path.exists(ssub_fn):
                ssub_d[pth] = open(ssub_fn, 'rb')
            else:
                raise ValueError('Runfolder {} does not contain a Sample-Submission-Form.xlsx'.format(pth))
    else:
        ssub_d[os.path.abspath(args.ssub.name)] = args.ssub
    args.ssub = ssub_d


    sample_dict = configmaker.merge_samples_with_submission_form(args.ssub,
                                                                 sample_dict,
                                                                 new_project_id=args.new_project_id,
                                                                 keep_batch=args.keep_batch)

    idat_dir = None
    if args.create_idat_dir:
        default_idat_dir = os.path.join('data', 'raw', 'idat')
        os.makedirs(default_idat_dir, exist_ok=True)
        s_ids = sample_dict.keys()
        if args.keep_batch:
            s_ids = set([s.split("_")[0] for s in s_ids])
        for sample_id in s_ids:
            for pid in project_dirs:
                r1_src, r2_src = match_idat(sample_id, pid, rel_path=False)
                r1_dst, r2_dst = match_idat(sample_id, pid, rel_path=True)
                if not any([r1_src, r1_dst, r2_src, r2_dst]):
                    continue
                for src, dst in zip(r1_src, r1_dst):
                    if src is not None:
                        dst = os.path.join(default_idat_dir, dst)
                        os.makedirs(os.path.dirname(dst), exist_ok=True)
                        os.symlink(src, dst)
                for src, dst in zip(r2_src, r2_dst):
                    if src is not None:
                        dst = os.path.join(default_idat_dir, dst)
                        os.makedirs(os.path.dirname(dst), exist_ok=True)
                        os.symlink(src, dst)
        idat_dir = default_fastq_dir


    config = configmaker.create_default_config(sample_dict, opts, args, fastq_dir=idat_dir)

    yaml.dump(config, args.output, default_flow_style=False, sort_keys=False)

    if args.create_project:
        if not os.path.exists("src/gcf-workflows"):
            os.makedirs('src', exist_ok=True)
            cmd = 'cd src && git clone {}'.format(GCF_WORKFLOWS_SRC)
            subprocess.check_call(cmd, shell=True)

        with open("src/gcf-workflows/libprep.config","r") as libprepconf_f:
            libconf = yaml.load(libprepconf_f, Loader=yaml.FullLoader)

        libkit = config['libprepkit']
        kitconf = libconf.get(libkit, None)
        if not kitconf:
            print("Libprepkit {} is not defined in libprep.config. Running with default settings.".format(libkit))
            workflow = "default"
        else:
            workflow = kitconf["workflow"]

        with open("Snakefile","w") as sn:
            sn.write(SNAKEFILE_TEMPLATE.format(workflow=workflow))


