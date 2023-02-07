#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
name: translate.py
date: 27.01.21
author: Snorre Sulheim

A script for translating metabolite IDs based on the metanetx data available from https://www.metanetx.org/mnxdoc/mnxref.html

This script is partially based and inspired by Daniel Machados translation class (metanetx.py) in the framed package
https://github.com/cdanielmachado/framed, doi: 10.5281/zenodo.1048261

"""

import pandas as pd
from pathlib import Path

import utils

DEFAULT_CHEM_XREF_PATH = utils.get_data_folder() / "metanetx" / "chem_xref.tsv"

class Translator(object):
    """ 
    Create draft metabolite ID translations based on the metanetx table chem_xref.tsv version 4.1
    """
    def __init__(self, chem_xref_path = None, remove_obsolete_IDs = True):
        self.set_chem_xref_path(chem_xref_path)
        self.parse_chem_xref(remove_obsolete_IDs)

    def set_chem_xref_path(self, chem_xref_path):
        if not chem_xref_path:
            self.chem_xref_path = DEFAULT_CHEM_XREF_PATH
        else:
            self.chem_xref_path = chem_xref_path

    def parse_chem_xref(self, remove_obsolete_IDs = True):
        """
        Read the chem_xref.tsv file and make a nicely formatted dataframe
        """
        df = pd.read_csv(self.chem_xref_path, sep = '\t', header = None, comment = '#', 
            names = ["external ID", "mnx ID", "Description"])

        # Drop obsolete IDs
        if remove_obsolete_IDs:
            # Many IDs are obsolete / redundant identifiers. These can be distinguished by their description 'secondary/obsolete/fantasy identifier'
            obsolete_rows = df["Description"] == "secondary/obsolete/fantasy identifier"
            df = df.loc[~obsolete_rows, :]

        # Drop the description column (currently not needed)
        df.drop('Description', axis=1, inplace = True)

        # Parse the mnx id (db:db_id) into seperate columns
        df[['db', 'db_id']] = df['external ID'].str.split(':', 1, expand=True)

        # Change a few db ids, e.g KeggC shoudl be kegg etc
        db_converter = {'keggC':'kegg', 'biggM':'bigg', 'seedM':'seed'}
        df["db"] = df["db"].replace(db_converter)
        self.met_df = df

    def translate_met(self, met_id, source_db, target_db, return_all = True, id_sorting = "alphabetically"):
        """
        Translate the metabolite ID from one db to another.
        Possible dbs:
            - kegg
            - bigg
            - slm
            - mnx
            - seed
            - chebi

        """
        source_db = source_db.lower()
        target_db = target_db.lower()
        db_list = ['kegg', 'bigg', 'slm', 'seed', 'chebi', 'mnx']

        for key, db in zip(["source DB", "target DB"], [source_db, target_db]):
            if not db in db_list:
                raise KeyError("The {0} is not in the list of possible databases: {1}".format(key, ", ".join(db_list)))


        # Find the MetaNetX IDs that correspond to the queried ID in the source DB
        if source_db == 'mnx':
            # match_rows = self.met_df['mnx_id'] == met_id
            # If the source ID is MetaNetX, we can basically skip the first search
            matched_mnx = set([met_id])
        else:
            match_rows = (self.met_df['db'] == source_db) & (self.met_df['db_id'] == met_id)
            # Collect the one (or possibly multiple() MetaNetX IDs that match the queried ID
            matched_mnx = set(self.met_df.loc[match_rows, 'mnx ID'])

        # Use the collected MetaNetX IDs to identify the correct ID in the target database
        if target_db == 'mnx':
            # If the target DB is MetaNetX we can skip this step
            target_IDs = matched_mnx
        else:
            target_IDs = []
            for mnx_id in matched_mnx:
                target_rows = (self.met_df['db'] == target_db) & (self.met_df['mnx ID'] == mnx_id)
                target_IDs += list(self.met_df.loc[target_rows, 'db_id'])
            target_IDs = list(set(target_IDs))

        # Sort mathed IDs
        if id_sorting == "alphabetically":
            target_IDs = sorted(target_IDs)
        elif id_sorting == "numerically":
            target_IDs = sorted(target_IDs, key = lambda x: re.findall(r"\d+", x)[0])
        else:
            pass

        if return_all:
            return target_IDs
        else:
            return target_IDs[0]





if __name__ == '__main__':
    T = Translator()

    T.translate_met('C00011', 'kegg', 'bigg')