#!/usr/bin/env python

import sys
import os
import re
import glob
import argparse


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
   matches = glob.glob(pth, 'SampleSheet.csv')
   return None

def inspect_dirs(args):
   project_dirs = []
   for pth in args.runfolders:
       pd = _match_project_dir(pth)
       project_dirs.append(pd)
   return project_dirs

def match_fastq(sample_name, project_dir):
   """Return fastq files matching a sample name.

   Returns paths relative to project directory
   """
   r1_fastq_files = glob.glob1(project_dir, '*' + sample_name + '*R1.fastq.gz', recursive=True)
   r2_fastq_files = glob.glob1(project_dir, '*' + sample_name + '*R2.fastq.gz', recursive=True)
   return None

def create_default_config(project_id=None):
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
   config['split']['skip'] = True
   config['split']['step'] = 'filter'
   config['split']['sscol'] = 'Sample_Project'
   return config


if __name__ == '__main__':

   parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument("-p", "--project-id" , help="Project ID", default="GCF-0000-000", type=is_valid_gcf_id)
   parser.add_argument("runfolders", nargs="+", help="Path(s) to flowcell dir(s)", action=FullPaths, type=is_dir)
   parser.add_argument("-s", "--sample-sheet", dest="samplesheet", type=argparse.FileType('r'), help="IEM Samplesheet")
   parser.add_argument("-o", "--output", default="config.yaml", help="Output config file", type=argparse.FileType('w'))
   parser.add_argument("--organism", default="homo_sapiens", help="Organism (if applicable to all samples)")
   parser.add_argument("--libkit", default="rna-trueseq", help="Library preparation kit. (if applicable for all samples)")
   parser.add_argument("-v", "--verbosity", action="count", default=0,
                       help="increase output verbosity")
   args = parser.parse_args()

   project_dirs = inspect_dirs(args)
   print(project_dirs)
