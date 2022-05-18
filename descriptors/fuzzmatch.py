"""
"""
import os
import sys
import requests
import copy
import pickle
import logging
import pprint

from thefuzz import fuzz, process
import pandas as pd

logger = logging.getLogger('GCF-configmaker')
logger.setLevel(10)

SERVER= "https://rest.ensembl.org"
ENS_ORG_DB = os.path.join(os.path.dirname(__file__), 'ens_org.pkl')
def fetch_ensembl_species(ext = "/info/species?", division='EnsemblVertebrates', db=None):
    kv_store = {}
    keep_keys = ['name', 'display_name', 'common_name', 'aliases']
    if division is not None:
        ext += 'division={}'.format(division)
    print(SERVER + ext)
    r = requests.get(SERVER + ext, headers={ "Content-Type" : "application/json"}) 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    species = decoded['species']
    for s in species:
        name = s.get('name')
        if name and '_' in name:
            genus_species = '_'.join(name.split('_')[:2])
        aliases = set(s.get('aliases', []))
        common_name = s.get('common_name')
        if common_name:
            aliases.add(common_name)
        strain = s.get('strain')
        if strain:
            aliases.add(strain)
        
        if None in aliases:
            aliases.remove(None)
        if genus_species != name:
            if strain is None:
                aliases.add(genus_species)
        kv_store[name] = list(aliases)
        
    print('found {} new species'.format(len(kv_store)))
    if db:
        db.update(kv_store)
        return db
    else:
        return kv_store

if os.path.exists(ENS_ORG_DB):
    with open(ENS_ORG_DB, 'rb') as fh:
        ORG_DB = pickle.load(fh)
else:
    ORG_DB = []
    for div in ['Vertebrates', 'Plants']:
        ORG_DB = fetch_ensembl_species(division='Ensembl'+ div, db=ORG_DB)
    with open(ENS_ORG_DB, 'wb') as fh:
        pickle.dump(ORG_DB, fh)


def fuzzmatch_organism(query, min_score=80, min_letters=4):
    """
    Fuzzy match organism name
    """
    if pd.isna(query) or query in ['N/A', 'NA', '', None]:
        return None
    best_score = 0
    best_match = {}
    logger.debug('Fuzzy matching org query: {}'.format(query))
    if query in ORG_DB.keys():
        logger.debug('genus_species exact match for {}.'.format(query))
        return query
    name, score = process.extractOne(query, ORG_DB.keys())
    if score >= min_score:
        logger.debug('fuzzymatch hit in genus_species key . {}, score: {}'.format(name, score))
        return name
    for org, choices in ORG_DB.items():
        score = 0
        filtered_choices = [i for i in choices if len(i) >= min_letters]
        if filtered_choices:
            filtered_choices.append(org)
            name, score = process.extractOne(query, filtered_choices)
        if score >= best_score:
            if score == best_score:
                best_match.append(org)
            else:
               best_match = [org] 
            best_score = score
    if len(best_match) > 1:
        logger.debug("several good hits ...")
        pprint.pprint(best_match, width=2)
        # ad hoc 1: remove strains and see if one remains
        no_strain_name = []
        for org in best_match:
            if org.count('_') == 1:
                no_strain_name.append(1)
            else:
                no_strain_name.append(0)
        if sum(no_strain_name) == 1:
            best_match = [best_match[no_strain_name.index(1)]]
            logger.debug("strain names removed. one hit remained: {}".format(best_match[0]))
        else:
            logger.debug('no strain names identification')

    if len(best_match) > 1:
        # ad hoc 2: keep best match with most aliases
        logger.debug("alias hitting: ")
        alias_hit, best_alias = False, ''
        best_hit, best_hit_score = None, 0
        for org in best_match:
            filtered_choices = [i for i in ORG_DB[org] if len(i) >= min_letters]
            matches = process.extractBests(org, filtered_choices, score_cutoff=min_score)
            n_hits = len(matches) if matches else 0
            if n_hits > 0:
                for name, score in matches:
                    if score > best_score:
                        best_match = [org]
                        best_score = score
                        alias_hit = True
                        best_alias = name
                if n_hits > best_hit_score:
                    best_hit_score = n_hits
                    best_hit = org

        if alias_hit:
            logger.debug("alias ({}) hit for org: {}".format(best_alias, best_match[0]))
        elif best_hit_score > 0:
            logger.debug('no aliases with higher score ')
            logger.debug('winning alias count ({}) for org: {}'.format(best_hit_score, org))
            best_match = [org]
        else:
            logger.debug('no hits in aliases ...')
                         
    if best_score >= min_score:
        if len(best_match)> 1:
            pprint.pprint(best_match, width=2)
        logger.debug('best match: {} ({})'.format(best_match[0], best_score))
        return best_match[0]
    else:
        logger.error("No acceptable matches by fuzzymatching organism (query=`{}`)".format(query))
        print("best matches: ")
        pprint.pprint(best_match, width=2)
        #logger.error(str(best_match[0]), best_score)
        return None

def fuzzmatch_reference(query, min_score=80):
    """
    Fuzzymatch common reference values
    """
    if pd.isnull(query) or query in ['N/A', 'NA', '', None]:
        return None
    keywords = ['control', 'ctrl', 'reference', 'wildtype', 'stable', 'zero', 'normal', 'healthy']
    name, score = process.extractOne(str(query).strip(), keywords)
    if score >= min_score:
        #print(name, score)
        logger.debug('identfied reference level {}'.format(query))
        return query
    return None
