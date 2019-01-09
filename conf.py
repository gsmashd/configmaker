#!/usr/bin/env python

import sys
import os
import re
import glob
import argparse
import pandas as pd
import yaml


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        values = [os.path.abspath(os.path.expanduser(v)) for v in values]
        setattr(namespace, self.dest, values)


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def is_valid_gcf_id(arg, patt='GCF-\d{4}-\d{3}'):
    m = re.match(patt, arg)
    if m:
        return m.group().strip()
    else:
        msg = "{0} is not a valid GCF number (format: GCF-YYYY-NNN)".format(arg)
        raise argparse.ArgumentTypeError(msg)

def _match_project_dir(pth):
    for fn in os.listdir(pth):
        if os.path.isdir(os.path.join(pth,fn)) and fn == args.project_id:
             return os.path.join(pth, fn)
    msg = "{0} is not present in run_folder: {1}".format(args.project_id, pth)
    raise argparse.ArgumentTypeError(msg)

def _match_samplesheet(pth):
    matches = glob.glob(os.path.join(pth, 'SampleSheet.csv'))
    return matches

def inspect_samplesheet(args):
    """
    if --samplesheet is set: Check that file exists and return it.
    else: check that runfolder(s) contain a SampleSheet.csv and return it (them).
    """
    if args.samplesheet is not None:
        if os.path.isfile(args.samplesheet):
            return args.samplesheet
        else:
            msg = "SampleSheet {} does not exist!".format(args.samplesheet)
            raise argparse.ArgumentTypeError(msg)
    else:
        samplesheets = []
        for pth in args.runfolders:
            ss = _match_samplesheet(pth)
            for s in ss:
                samplesheets.append(s)
        if len(samplesheets) == 0:
            msg = "Cannot find SampleSheet.csv in {}".format(', '.join(args.runfolders))
            raise RuntimeError(msg)
        return samplesheets

def get_project_samples_from_samplesheet(args):
    """
    Return a dataframe containing project samples
    """
    ss = inspect_samplesheet(args)
    df_list = []
    with open(ss[0],'r') as s:
        while True:
            line = s.readline()
            if not line:
                msg = 'No [data]-section in samplesheet {}'.format(s.name)
                raise RuntimeError(msg)
            if line.startswith('[Data]'):
                df_list.append(pd.read_csv(s))
                break
    for i,ss_df in enumerate(df_list):
        df_list[i] = ss_df[ss_df.Sample_Project == args.project_id]
    samples = {}
    for df in df_list:
        for sample_d in df.to_dict('records'):
            sample = samples.get(sample_d['Sample_ID'],None):
            if sample:
                if sample == sample_d:
                    continue
                else:
                    msg = 'The samplesheets contains conflicting records for Sample_ID {}'.format(sample_d['Sample_ID'])
                    raise RuntimeError(msg)
            samples[sample_d['Sample_ID']] = sample_d

    return pd.DataFrame.from_dict(samples,orient='index')

def inspect_dirs(args):
    project_dirs = []
    for pth in args.runfolders:
        pid = _match_project_dir(pth)
        project_dirs.append(pid)
    return project_dirs

def match_fastq(sample_name, project_dir):
    """Return fastq files matching a sample name.

    Returns paths relative to project directory
    """
    r1_fastq_files = glob.glob(os.path.join(project_dir, sample_name + '*R1.fastq.gz'), recursive=True)
    r2_fastq_files = glob.glob(os.path.join(project_dir, sample_name + '*R2.fastq.gz'), recursive=True)

    r1 = [os.path.relpath(x,os.path.dirname(os.path.dirname(project_dir))) for x in r1_fastq_files]
    r2 = [os.path.relpath(x,os.path.dirname(os.path.dirname(project_dir))) for x in r2_fastq_files]
    return r1, r2

def find_samples(df, project_dirs):
    sample_dict = {}
    for index, row in df.iterrows():
        s_r1 = []
        s_r2 = []
        for p_pth in project_dirs:
            r1, r2 = match_fastq(row.Sample_ID, p_pth)
            s_r1.extend(r1)
            s_r2.extend(r2)
        pe = 0 if len(s_r2) == 0 else 1
        sample_dict[row.Sample_ID] = {
                'R1': ','.join(s_r1),
                'R2': ','.join(s_r2),
                'paired_end': pe,
                'Sample_ID': row.Sample_ID,
            }

    return sample_dict

def create_default_config(sample_dict,project_id=None):
    config = {}
    if project_id:
         config['project_id'] = project_id
    config['ext_dir'] = 'data/ext'
    config['interim_dir'] = 'data/tmp'
    config['processed_dir'] = 'data/processed'
    config['merge'] = {}
    config['merge']['skip'] = False
    config['merge']['step'] = 'quant'
    config['merge']['sscol'] = 'Sample_ID'
    config['split'] = {}
    config['split']['skip'] = True
    config['split']['step'] = 'filter'
    config['split']['sscol'] = 'Sample_Project'
    config['samples'] = sample_dict
    return config


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--project-id" , help="Project ID", default="GCF-0000-000", type=is_valid_gcf_id)
    parser.add_argument("runfolders", nargs="+", help="Path(s) to flowcell dir(s)", action=FullPaths, type=is_dir)
    parser.add_argument("-s", "--sample-sheet", dest="samplesheet", type=argparse.FileType('r'), help="IEM Samplesheet")
    parser.add_argument("-o", "--output", default="config.yaml", help="Output config file", type=argparse.FileType('w'))
    parser.add_argument("--sample-submission-form", dest="ssub", type=argparse.FileType('r'), help="GCF Sample Submission Form")
    parser.add_argument("--organism", default="homo_sapiens", help="Organism (if applicable to all samples)")
    parser.add_argument("--libkit", default="rna-trueseq", help="Library preparation kit. (if applicable for all samples)")
    parser.add_argument("-v", "--verbosity", action="count", default=0,
                              help="increase output verbosity")
    args = parser.parse_args()

    project_dirs = inspect_dirs(args)
    s_df = get_project_samples_from_samplesheet(args)
    sample_dict = find_samples(s_df,project_dirs)
    config =  create_default_config(sample_dict,args.project_id)
    yaml.dump(config,args.output,default_flow_style=False)

