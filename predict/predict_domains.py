import os, sys
import pickle, csv, json
from collections import *
from itertools import *
from tqdm import tqdm as tqdm

import numpy as np

import utils

# Models are specified using the named tuple (in utils/load_models.py):
# ModelsSpecificationFormat = namedtuple("ModelsSpecificationFormat", ["models_format", "models"]) 

def predict_domain_peptide_interactions(domain_file, peptide_file, models, output_file, progressbar=False):
    """
    Predicts domain peptide interactions between the domain and peptide metadata files.
    """
    
    domain_metadata = list(csv.reader(open(domain_file, 'r')))
    peptide_metadata = list(csv.reader(open(peptide_file, 'r')))

    with open(output_file, 'w+') as ofile:
        writer = csv.writer(ofile, delimiter=',')
        
        total = len(domain_metadata) * len(peptide_metadata)
        for (did, dseq, dtype), (pid, pseq, ptype, mtype) in tqdm(product(domain_metadata, peptide_metadata), disable=(not progressbar), desc="Pairwise Domain-peptide iteration", total=total):
            if dtype in models.models and ptype in models.models[dtype]:
                v = models.models[dtype][ptype](dseq, pseq)
                writer.writerow([dtype, did, dseq, ptype, pid, pseq, v])

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("domains_metadata", type=str)
    parser.add_argument("peptides_metadata", type=str)
    parser.add_argument("-o", "--output_likelihoods", type=str, default="outputs/domain_peptide_predictions.csv")
    parser.add_argument("-p", "--progressbar", action='store_true', default=False) 

    model_specification_group = parser.add_argument_group(title="Trained models specification")
    model_specification_group.add_argument("-m", "--models", type=str, default="models/hsm_pretrained/", 
            help="Defines a directory containing ." )
    model_specification_group.add_argument("--model_format", type=str, default="models/hsm_pretrained/model_formats.csv",
            help="Defines the parameters for each model. Default: models/hsm_pretrained/model_formats.csv")
    model_specification_group.add_argument("--amino_acids_order", type=str, default="models/amino_acid_ordering.txt",
            help="Define a different group (and order) of amino-acids. Needed to add a 'chemistry' type not currently used " +
            "like phospho-serine/threonine. Default: models/amino_acid_ordering.txt.")
    
    opts = parser.parse_args()
    
    models_specification = utils.load_models.load_models_from_dir(opts.models, opts.model_format, opts.amino_acids_order)
    predict_domain_peptide_interactions(opts.domains_metadata, opts.peptides_metadata, models_specification, opts.output_likelihoods, progressbar=opts.progressbar)
