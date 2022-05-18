#!/usr/bin/env python
"""descriptors are key-value pairs to describe properties of a columns in a sample table

descritpors are stored in a dictionary of dictionary with main main dict keyed on sample_id and subdict is key:value

supported descriptors:

dtype: accepted items: category, numerical, unique (nominal, ordinal, continous, discrete)
ref_level: string value (must have dtype category and be a valid category)
min: number (must have numerical dtype)
max: number (must have numerical dtype)
mean_var: [mean, var] theoretical/emprical mean and variance
title:
description:
suffix:
format:
modify:
shared_key:
scale:
placement:
model_pkl:


Descriptors are coded in headers within brackets by comma separated key=value pairs
examples:
Sample_Group[dtype=category, ref_level=control]
RIN [dtype=numerical, min=0, max=10]

Descriptors may also be added by a yaml

"""
DTYPES = ['categorical', 'numerical', 'string', 'unique', 'enum', 'boolean', 'constant', 'nucleotides']

import os
import re
import logging
import string
import oyaml as yaml

import pandas as pd

from descriptors import fuzzmatch

logger = logging.getLogger('GCF-configmaker')

DEFAULT_DESCRIPTORS = {}
desc_fn = os.path.join(os.path.dirname(__file__), "default_descriptors.yaml")
if os.path.exists(desc_fn):
    with open(desc_fn) as fh:
       DEFAULT_DESCRIPTORS  = yaml.load(fh, Loader=yaml.FullLoader)


def findall_header_descriptors(df, strip_descriptors=True, mapper=None):
    """
    Extract decriptor tags in dataframe column names.

    Parameters
    ----------
    df : dataframe
        sample table
    strip_descriptors : boolean
        Remove descriptor tags from datfarme column names
    mapper : dict-like or function
        Dict-like or functions transformations to apply to column names

    Returns
    -------
    dict
        Diconary of dictionary. Descriptor key-value pairs for eacn column name 
    
    
    """
    patt = re.compile('\[(.+?)\]')
    desc = {}
    for e in df.columns:
        e = str(e)
        name = e.split(r'[')[0].strip()
        m = patt.search(e)
        if m:
            kv = {}
            match_str = patt.findall(e)
            for s in match_str:
                match_list = [str(i).strip() for i in s.split(',')]
                for el in match_list:
                    k, v = el.split('=')
                    kv[str(k).strip()] = str(v).strip()
        else:
            kv = {}
        if mapper is not None:
            name = mapper(name)
        desc[name] = kv
    if strip_descriptors:
        df.columns =  [i.split(r'[')[0].strip() for i in df.columns]
    return desc


def add_default_descriptors(df, desc):
    """
    Add default descriptor key-value pairs
    
    Params
    ------
    df : dataframe
        sample table
    desc : dict
        Dictionary of descriptors

    Returns
    -------
    dict
        Updated dicionary of descriptors
        
    """
    unknown = DEFAULT_DESCRIPTORS['unknown']
    for col_name in df.columns:
        if not col_name in desc:
            desc[col_name] = {}
        for k, v in DEFAULT_DESCRIPTORS.get(col_name, unknown).items():
            if k not in desc[col_name]:
                logger.debug('default descriptor ({}:{}) to {}'.format(k, v, col_name))
                desc[col_name][k] = v
        if col_name not in desc:
            logger.debug('mssing default descriptor for {}'.format(col_name))
    return desc


def order_columns_by_descriptors(df, desc):
    logger.debug('reordering columns by descriptors ...')
    placement = []
    for col in df.columns:
        p = desc.get(col, {}).get('placement', 1000)
        placement.append(p)
    order = [i[0] for i in sorted(enumerate(placement), key=lambda x:x[1])]
    df = df.iloc[:,order]
    return df



def _is_part_numerical(col):
    n_numeric = col.str.isnumeric().sum()
    if (n_numeric > 0) and (n_numeric < len(col)):
        return True
    return False
    

def write_descriptors(desc, fn):
    """
    Serialize descriptors to YAML format
    """
    with open(fn, "w") as fh:
        yaml.dump(desc, fh, Dumper=yaml.Dumper)

DTYPES = ['categorical', 'numerical', 'string', 'unique', 'enum', 'boolean', 'constant', 'nucleotide_sequence', 'date', 'path']
SUBTYPES = ['unique', 'enum', '']



def _infer_dtype_string(x, subtype=None, col_name=None):
    """
    accepted string values are alphanumeric + underscore/whitespace/comma - punctuation
    norwegian letters are translated
    """
    x = str(x).strip()
    if subtype == 'no_conversion':
        return x
    elif subtype == 'gcf_number':
        m = re.match('GCF-\d{4}-\d{3}', x)
        if m:
            return m.group().strip()
        else:
            msg = "{0} is not a valid GCF number (format: GCF-YYYY-NNN)".format(x)
            logger.error(msg)
    elif subtype == 'nucleotide':
        if x not in ['C', 'T', 'A', 'G']:
            logger.error('{}: illegal nucleotide letter {}'.format(name, x))
        
    x = ''.join(filter(lambda x: x.replace('_', '').replace(' ', '').isalnum, x))
    trans = str.maketrans('', '', string.punctuation.replace('_', '').replace(',', ''))
    x = x.translate(trans)
    x = x.replace('ø', 'oe').replace('Ø', 'OE')
    x = x.replace('å', 'aa').replace('Å', 'AA')
    x = x.replace('æ', 'ae').replace('Æ', 'AE')
    
    return x


def _tryhard_numeric(x):
    """
    """
    KNOWN_PREFIXES = ['>', '<', '>=', '<=', 'plate', '~', '-']
    RANGE_PATT = '([\d\.]+)\s?-\s?([\d\.]+)'
    DELIMTERS = ['|', '/', '&']
    units = u'(?:mg|µg|ng|ml|µl|nl|mM|µM|nM|ng/µl)'
    PATT = '(^\D*)([\d\.]+)(\D*$)'

    logger.debug('numerical tryhard conversion ...')
    x = str(x).strip()
    r = re.search(RANGE_PATT, x)
    if r:
        logger.debug('identified a range pattern {}'.format(x))
        vals = r.groups()
        mean = sum(map(float, vals))/float(len(vals))
        logger.debug('identified a range pattern {}: mean={}'.format(x, mean))
        return mean
    for d in DELIMTERS:
        els = x.split(d)
        if len(els) > 1:
            if all([re.match('[.\d]+', i) for i in els]):
                logger.debug('identified a delimited pattern {}'.format(x))
                return '|'.join(els)    
    m = re.search(PATT, x)
    if m:
        prefix, num, postfix = [i.strip() for i in m.groups()]
        logger.debug('numerical tryhard conversion pattern prefix={},num={},postfix={}'.format(prefix, str(num), postfix))
        return float(num)
    logger.debug('no th conversion.')
    return x
    
            


def _infer_dtype_numerical(x, tryhard=False, subtype=None, force_int=True, col_name=None):
    """
    numerical values are decimal or integers
    """
    if pd.isnull(x) or not x:
        return x
    xx = str(x).replace(',', '.')
    negative_sign = True if xx.startswith('-') else False
    if negative_sign:
        xx = xx[1:]
    out = None
    try:
        if xx.replace(' ', '').replace('.', '').isdigit():
            if subtype == 'float':
                out = float(xx)
            elif subtype == 'int':
                #if int(float(x)) == float(x):
                out =  int(float(xx)) # downcast
                #else:
                #    raise ValueError
            elif '.' in xx:
                out = float(xx)
            else:
                out = int(xx)
        else:
            raise ValueError
    except:
        if tryhard:
            xx = _tryhard_numeric(xx)
            if xx:
                out = _infer_dtype_numerical(xx, subtype=subtype, tryhard=False, col_name=col_name)
    
    if out: # conversion ok
        if negative_sign:
            out = -out
        return out
    else:
        subtype = subtype or 'numerical'
        name = col_name or ''
        logger.error('{}: failed to sanitize to {} `{}`'.format(name, subtype, str(x)))
        return x
        
def _infer_dtype_categorical(x, subtype=None, enum=None, col_name=None):
    if subtype == 'int':
        out = _infer_dtype_numerical(x, tryhard=True, subtype='int', col_name=col_name)
    elif subtype == 'enum':
        # fuzzmatch close ?
        if enum:
            assert(x in enum)
        else:
            out = x
    else:
        out = _infer_dtype_string(x, subtype='no_conversion', col_name=col_name)
    return out




def infer_by_descriptor(df, desc):
    df = df.convert_dtypes()
    for col_name in df.columns:
        col = df[col_name]
        dtype = desc.get(col_name, {}).get('dtype')
        subtype = desc.get(col_name, {}).get('subtype')
        if not col_name in desc:
            desc[col_name] = {}
        if dtype == 'numerical':
            col = col.apply(_infer_dtype_numerical, tryhard=True, subtype=subtype, col_name=col_name)
            col = pd.to_numeric(col, errors='coerce')
            col_min = desc.get(col_name, {}).get('min', col.min())
            if col_min > col.min():
                logger.warning('decriptor min value > data min value {}>{}'.format(col_min, col.min()))
            desc[col_name]['min'] = float(col_min)
            col_max = desc.get(col_name, {}).get('max', col.max())
            if col_max < col.max():
                logger.warning('decriptor max value < data max value {}<{}'.format(col_max, col.max()))
            desc[col_name]['max'] = float(col_max)
            if subtype is None:
                if col.dtype.kind == 'i':
                    subtype = 'int'
                    col = pd.to_numeric(col, downcast='integer', errors='coerce')
                else:
                    subtype = 'float'
                    col = pd.to_numeric(col, downcast='float', errors='coerce')
        elif dtype == 'string':
            col = col.apply(_infer_dtype_string, subtype=subtype, col_name=col_name)
            if subtype == 'unique':
                assert(len(set(col)) == len(col))
            elif subtype == 'organism':
                org_map = {}
                for org in set(col):
                    fm = fuzzmatch.fuzzmatch_organism(org)
                    if fm:
                        org_map[org] = fm
                for i, c in enumerate(col):
                    col[i] = org_map.get(c) or pd.NA
            elif subtype is None:
                if len(set(col)) == 1:
                    subtype = 'constant'
                elif len(set(col)) == len(col):
                    subtype = 'unique'
            col = col.astype('string')
            

        elif dtype == 'categorical':
            col = col.apply(_infer_dtype_categorical, subtype=subtype, col_name=col_name)
            col = col.fillna('')
            enum = desc.get(col_name, {}).get('enum', list(set(col)))
            ordered = desc.get(col_name, {}).get('ordered', False)
            
            if ordered:
                if '' in enum:
                    enum.remove('')
                if subtype == 'int':
                    enum = sorted(enum, key=lambda x: x if str(x).isdigit() else 99999999)
                else:
                    enum = sorted(enum)
            if '' in col and '' not in enum:
                enum.append('')
            
            cat_type = pd.api.types.CategoricalDtype(categories=enum, ordered=ordered)
            ref_level = desc.get(col_name, {}).get('ref_level')
            if ref_level is None:
                if len(cat_type.categories) > 1:
                    for v in cat_type.categories:
                        query = str(v).strip()
                        logger.debug('{}: fuzzymatching {} as reference level ...'.format(col_name, query))
                        ref_level = fuzzmatch.fuzzmatch_reference(query)
            if ref_level and ref_level in cat_type.categories:
                enum = list(cat_type.categories)
                enum.remove(ref_level)
                enum.insert(0, ref_level)
                cat_type = pd.api.types.CategoricalDtype(categories=enum, ordered=ordered)
                desc[col_name]['ref_level'] = ref_level
            col = col.astype(cat_type)
            if len(cat_type.categories) == 1:
                subtype = 'constant'
            else:
                if col_name in ['Sample_Biosource', 'Sample_Group']:
                    subtype = subtype or 'experimental_value'
        elif dtype == 'boolean':
            col = col.astype(dtype)
        elif dtype == 'nucleotide_sequence':
            col = col.astype('string').str.upper()
        elif dtype in DTYPES or dtype is None:
            pass
        else:
            raise ValueError('unknown dtype {}'.format(dtype))
        if subtype:
            desc[col_name]['subtype'] = subtype
        df[col_name] = col
    return order_columns_by_descriptors(df, desc), desc

