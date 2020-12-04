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



SEQUENCERS = {
    'NB501038' : 'NextSeq 500',
    'SN7001334' : 'HiSeq 2500',
    'K00251' : 'HiSeq 4000',
    'M02675' : 'MiSeq NTNU',
    'M03942' : 'MiSeq StOlav',
    'M05617' : 'MiSeq SINTEF'
}

PIPELINE_MAP = {
    'Lexogen SENSE Total RNA-Seq Library Prep Kit (w/RiboCop rRNA Depletion Kit V1.2)': 'rna-seq',
    'Lexogen SENSE mRNA-Seq Library Prep Kit V2': 'rna-seq',
    'Illumina TruSeq Stranded Total RNA Library Prep (Human/Mouse/Rat)': 'rna-seq',
    'Illumina TruSeq Stranded Total RNA Library Prep (Globin)': 'rna-seq',
    'Illumina TruSeq Stranded mRNA Library Prep': 'rna-seq',
    'QIAseq 16S ITS Region Panels': 'microbiome',
    '16S Metagenomic Sequencing Library Prep': 'microbiome',
    'ITS Low Input GCF Custom': 'microbiome',
    '10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3': 'single-cell',
    'Bioo Scientific NEXTflex Small RNA-Seq Kit v3': 'small-rna'
}

REPO_MAP = {
    'rna-seq': 'https://github.com/gcfntnu/rna-seq.git',
    'single-cell': 'https://github.com/gcfntnu/single-cell.git',
    'microbiome': 'https://github.com/gcfntnu/microbiome.git',
    'small-rna': 'https://github.com/gcfntnu/small-rna.git'
}

GCFDB_SRC = "https://github.com/gcfntnu/gcfdb.git"

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
    if arg is None:
        return True
    m = re.match(patt, arg)
    if m:
        return m.group().strip()
    else:
        msg = "{0} is not a valid GCF number (format: GCF-YYYY-NNN)".format(arg)
        raise argparse.ArgumentTypeError(msg)

def _match_project_dir(pth, project_id=None):

    if project_id:
        for fn in os.listdir(pth):
            if os.path.isdir(os.path.join(pth, fn)) and fn in project_id:
                return os.path.join(pth, fn), fn
        msg = "{0} is not present in run_folder: {1}".format(project_id, pth)
        raise ValueError(msg)
    else:
        project_dir = None
        for fn in os.listdir(pth):
            if os.path.isdir(os.path.join(pth, fn)) and re.match('GCF-\d{4}-\d{3}', fn):
                if project_dir is not None:
                    raise ValueError('runfolders contain more than one project folders. Use `--project-id` option to choose one.')
                project_dir = os.path.join(pth, fn)
                project_id = fn
        if project_dir:
            return project_dir, project_id
        raise ValueError('failed to identify any valid projects in runfolder: {}'.format(pth))

def _match_samplesheet(pth):
    matches = glob.glob(os.path.join(pth, '*SampleSheet*.csv'))
    return matches

def inspect_samplesheet(samplesheet, runfolders):
    """
    if --samplesheet is set: Check that file exists and return it.
    else: check that runfolder(s) contain a SampleSheet.csv and return it (them).
    """
    if samplesheet is not None:
        return [samplesheet.name]
    else:
        samplesheets = []
        for pth in runfolders:
            ss = _match_samplesheet(pth)
            for s in ss:
                samplesheets.append(s)
        if len(samplesheets) == 0:
            msg = "Cannot find SampleSheet.csv in {}".format(', '.join(runfolders))
            raise RuntimeError(msg)
        return samplesheets

def get_data_from_samplesheet(fh):
    custom_opts = False
    opts_d = {}
    while True:
        line = fh.readline()
        if not line:
            msg = 'No [data]-section in samplesheet {}'.format(s.name)
            raise RuntimeError(msg)
        if line.startswith('[Data]'):
            return pd.read_csv(fh), opts_d
        elif line.startswith('[CustomOptions]'):
            custom_opts = True
            continue
        elif custom_opts:
            key, val = [i.rstrip() for i in line.split(',')]
            if val.lower() == 'true':
                val = True
            opts_d[key] = val

def get_project_samples_from_samplesheet(samplesheet, runfolders, project_id):
    """
    Return a dataframe containing project samples
    """
    ss = inspect_samplesheet(samplesheet, runfolders)
    df_list = []
    for sheet in ss:
        with open(sheet, 'r') as s:
            data, opts = get_data_from_samplesheet(s)
            df_list.append(data)
    df = pd.concat(df_list)
    if project_id:
        df = df[df.Sample_Project.isin(project_id)]
    df['Sample_ID'] = df['Sample_ID'].astype(str)
    df = df[['Sample_ID']]
    df = df.drop_duplicates(['Sample_ID'])
    return df, opts

def inspect_dirs(runfolders, project_id=None):
    project_dirs = []
    project_ids = set()
    for pth in runfolders:
        pdir, pid = _match_project_dir(pth, project_id)
        project_dirs.append(pdir)
        project_ids.add(pid)
    if len(project_ids) == 0:
        raise ValueError('runfolders does not contain any of the specified projects.')
    return project_dirs, list(project_ids)

def match_fastq(sample_name, project_dir, rel_path=True):
    """Return fastq files matching a sample name.

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
        warnings.warn(warn_msg)
        return None, None
    r1_fastq_files = sorted(r1_fastq_files)
    r2_fastq_files = sorted(r2_fastq_files)
    if rel_path:
        r1_fastq_files = [os.path.relpath(x, os.path.dirname(os.path.dirname(project_dir))) for x in r1_fastq_files]
        r2_fastq_files = [os.path.relpath(x, os.path.dirname(os.path.dirname(project_dir))) for x in r2_fastq_files]

    return r1_fastq_files, r2_fastq_files

def find_samples(df, project_dirs):
    sample_dict = {}
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
            warnings.warn(warn_str)
        else:
            pe = 0 if len(s_r2) == 0 else 1
            sample_dict[str(row.Sample_ID)] = {
                'R1': ','.join(s_r1),
                'R2': ','.join(s_r2),
                'paired_end': pe,
                'Sample_ID': row.Sample_ID,
            }
    return sample_dict


def find_samples_batch(df, project_dirs):
    sample_dict = {}
    for index, row in df.iterrows():
        for p_pth in project_dirs:
            r1, r2 = match_fastq(row.Sample_ID, p_pth)
            if (not r1) and (not r2):
                warn_str = 'sample {} not found in {}'.format(row.Sample_ID, p_pth)
                warnings.warn(warn_str)
            else:
                pe = 0 if not r2 else 1
                r2 = [] if not r2 else r2
                sample_dict[str(row.Sample_ID) + "_" + p_pth.split("/")[-2].split("_")[-1]] = {
                    'R1': ','.join(r1),
                    'R2': ','.join(r2),
                    'paired_end': pe,
                    'Sample_ID': str(row.Sample_ID) + "_" + p_pth.split("/")[-2].split("_")[-1],
                    'Src_Sample_ID': row.Sample_ID,
                }
    return sample_dict

def sample_submission_form_parser(ssub_path, keep_batch=None):
    customer_column_map = {
        'Unique Sample ID': 'Sample_ID',
        'External ID (optional reference sample ID)': 'External_ID',
        'Sample Group (conditions to be compared)': 'Sample_Group',
        'Comments (optional info that does not fit in other columns)': 'Customer_Comment',
        'Sample biosource (examples: celltype/tissue/FFPE)': 'Sample_Biosource',
        'Project ID': 'Project_ID',
        'Sample type (e.g RNA or DNA or library)': 'Sample_Type',
        'Index2_p5 (If dual indexed libraries are submitted, indicate what index sequence is used as p5)': 'Index',
        'Index1_p7 (If libraries are submitted, indicate what index sequence is used as P7)': 'Index2',
        'Plate location (if samples delivered in 96 well plates)': 'Plate',
        'Sample Buffer': 'Sample_Buffer',
        'Volume (ul)': 'Volume',
        'Quantification Method': 'Quantification_Method',
        'Concentration (ng/ul)': 'Concentration',
        '260/280 ratio': '260/280',
        '260/230 ratio': '260/230',
        }
    lab_column_map = {
            'Concentration (ng/ul)': 'Concentration',
            '260/280 ratio': '260/280',
            '260/230 ratio': '260/230',
            'Comment': 'Lab_Comment'
        }

    customer = pd.read_excel(ssub_path, sheet_name=0, skiprows=14)
    customer.rename(columns=customer_column_map, inplace=True, errors='ignore')
    remove_cols = ['Concentration', 'Index', 'Index2', 'Sample_Type', 'Plate', 'Sample_Buffer', 'Volume', 'Quantification_Method', 'Concentration', '260/280', '260/230']
    customer = customer.drop(remove_cols, axis=1, errors='ignore')

    lab = pd.read_excel(ssub_path, sheet_name=2)
    lab.rename(columns=lab_column_map, inplace=True, errors='ignore')
    lab.drop(['Sample_Name','Project ID','KIT'], inplace=True, errors='ignore')

    if keep_batch:
        customer['Sample_ID'] = customer['Sample_ID'].astype(str) + "_" + os.path.dirname(ssub_path).split("_")[-1]
        customer['Flowcell_Name'] = [os.path.basename(pth)]*len(customer)
        customer['Flowcell_ID'] = [os.path.basename(pth).split("_")[-1]]*len(customer)
        lab['Sample_ID'] = lab['Sample_ID'].astype(str) + "_" + os.path.dirname(ssub_path).split("_")[-1]
    if not lab.empty:
        merge_ssub = pd.merge(customer, lab, on='Sample_ID', how='inner')
    else:
        merge_ssub = customer
    merge_ssub["Sample_ID"] = merge_ssub["Sample_ID"].astype(str)
    merge_ssub.index = merge_ssub["Sample_ID"]
    merge_ssub_d = merge_ssub.to_dict(orient='index')
    return merge_ssub_d

def merge_samples_with_submission_form(ssub, sample_dict, new_project_id=None, keep_batch=None):
    merge_d = dict()
    for pth in ssub.keys():
        merge_ssub_d = sample_submission_form_parser(ssub[pth].name, keep_batch=keep_batch)
        merge_d[pth] = merge_ssub_d

    merge = dict()
    for pth, s_d in merge_d.items():
        intersection = set(merge.keys()).intersection(set(s_d.keys()))
        if len(intersection) > 0:
            logger.warning("WARNING: Sampleinfo on {} are updated with values from {}/Sample-Submission-Form.xlsx. Specify a custom sample submission form with --sample-submission-form to force values.".format(', '.join(list(intersection)), pth))
        merge.update(s_d)
    merge = pd.DataFrame.from_dict(merge, orient='index')
    check_existence_of_samples(sample_dict.keys(), merge)
    sample_df = pd.DataFrame.from_dict(sample_dict, orient='index')
    sample_df = sample_df.merge(merge, on='Sample_ID', how='left')
    sample_df.reset_index()
    sample_df.index = sample_df['Sample_ID']
    sample_df.fillna('', inplace=True)
    if new_project_id:
        sample_df.rename(columns={'Project_ID': 'Src_Project_ID'}, inplace=True)
        sample_df.insert(loc=0, column="Project_ID", value=[new_project_id]*len(sample_df))
    s_dict = sample_df.to_dict(orient='index')
    return s_dict

def check_existence_of_samples(samples, df):
    diff = set(samples) - set(df['Sample_ID'].astype(str))
    if diff:
        logger.warning("WARNING: Samples {} are contained in SampleSheet, but not in sample submission form. Sample info from these samples will have empty values.".format(', '.join(list(diff))))
    diff = set(df['Sample_ID'].astype(str)) - set(samples)
    if diff:
        logger.warning("WARNING: Samples {} are contained in sample submission form, but not in SampleSheet. Sample info from these samples are omitted.".format(', '.join(list(diff))))
    return None

def find_read_geometry(runfolders):
    all = set()
    for fn in runfolders:
        stats_fn = os.path.join(fn, 'Stats', 'Stats.json')
        read_geometry = []
        with open(stats_fn) as fh:
            S = json.load(fh)
        for read in S['ReadInfosForLanes'][0]['ReadInfos']:
            if not read['IsIndexedRead']:
                read_geometry.append(read['NumCycles'])
        all.add(':'.join(map(str, read_geometry)))
    if len(all) > 1:
        raise ValueError('Read geometry mismatch between runfolders. Check Stats.json!')
    return read_geometry

def find_machine(runfolders):
    matches = set()
    for pth in runfolders:
        machine_code = os.path.basename(pth).split('_')[1]
        machine = SEQUENCERS.get(machine_code, '')
        matches.add(machine)
    if len(matches) > 1:
        logger.warning('Multiple sequencing machines identified!')
    return '|'.join(list(matches))

def create_default_config(sample_dict, opts, args, fastq_dir=None):
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

    config['samples'] = sample_dict

    return config


if __name__ == '__main__':

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
    parser.add_argument("--keep-batch", action='store_true', help="Sample names will be made unique for each batch.")

    args = parser.parse_args()

    project_dirs, args.project_id = inspect_dirs(args.runfolders, args.project_id)
    s_df, opts = get_project_samples_from_samplesheet(args.samplesheet, args.runfolders, args.project_id)
    if args.keep_batch:
        sample_dict = find_samples_batch(s_df, project_dirs)
    else:
        sample_dict = find_samples(s_df, project_dirs)

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


    sample_dict = merge_samples_with_submission_form(args.ssub,
                                                     sample_dict,
                                                     new_project_id=args.new_project_id,
                                                     keep_batch=args.keep_batch)

    fastq_dir = None
    if args.create_fastq_dir:
        default_fastq_dir = os.path.join('data', 'raw', 'fastq')
        os.makedirs(default_fastq_dir, exist_ok=True)
        s_ids = sample_dict.keys()
        if args.keep_batch:
            s_ids = set([s.split("_")[0] for s in s_ids])
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
        fastq_dir = default_fastq_dir


    config = create_default_config(sample_dict, opts, args, fastq_dir=fastq_dir)

    yaml.dump(config, args.output, default_flow_style=False, sort_keys=False)

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
    with open(os.path.join(dirname, ".configmaker.log"),"w") as conflog:
        print("Summary:")
        for nf, ns in count.items():
            line = "{} sample{} found in {} flowcell{}".format(ns, "s" if ns > 1 else "", nf, "s" if nf > 1 else "")
            print(line)
            conflog.write(line + "\n")
        conflog.write("Sample summary:\n")
        for s, f in summary.items():
            conflog.write("Sample {} found in: {}\n".format(s, ', '.join(f)))
    print("Full sample summary log written to {}".format(os.path.join(dirname, ".configmaker.log")))


    if args.create_project:
        pipeline = PIPELINE_MAP.get(config.get('libprepkit', None), None)
        if pipeline:
            if not os.path.exists("src/{}".format(pipeline)):
                os.makedirs('src', exist_ok=True)
                cmd = 'cd src && git clone {}'.format(REPO_MAP.get(pipeline, None))
                subprocess.check_call(cmd, shell=True)

            if not (os.path.exists("src/gcfdb") or os.environ.get("GCF_DB")):
                cmd = 'cd src && git clone {}'.format(GCFDB_SRC)
                subprocess.check_call(cmd, shell=True)

            cmd = 'wget -O Snakefile https://gcf-winecellar.medisin.ntnu.no/snakefiles/Snakefile-{}'.format(pipeline)
            subprocess.check_call(cmd, shell=True)
        else:
            raise ValueError('Libprepkit {} is not associated with any Snakemake pipelines'.format(config.get('libprepkit')))

