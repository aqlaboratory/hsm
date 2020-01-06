import os, sys
import pickle, csv, json
from collections import *
from itertools import *
import tqdm

import numpy as np

import utils

# Models are specified using the named tuple (in utils/load_models.py):
# ModelsSpecificationFormat = namedtuple("ModelsSpecificationFormat", ["models_format", "models"]) 

MetadataTuple = namedtuple("MetadataTuple", ["pairs", "domain_metadata", "peptide_metadata"])
def load_metadata(domain_metadata_fpath, peptide_metadata_fpath, pretrained_models, ppi_pairs_fpath=None):
    """
    Loads input metadata from the pre-trained models.
    """
    
    domain_mdata = pickle.load(open(domain_metadata_fpath, 'rb'))
    peptide_mdata = pickle.load(open(peptide_metadata_fpath, 'rb'))
   
    if ppi_pairs_fpath is None:
        pairs = set(tuple(sorted(t)) for t in product(set(domain_mdata.keys()), set(peptide_mdata.keys())))
    else:
        pairs = set(tuple(sorted(r)) for r in csv.reader(open(ppi_pairs_fpath, 'r')))
    
    pairs = [t for t in pairs if utils.ppi_prediction.is_valid(*t, domain_mdata, peptide_mdata, pretrained_models)]

    return MetadataTuple(pairs, domain_mdata, peptide_mdata)

def predict_interactions(metadata, models, output_ppis_fname, output_dpis_fname):
    """
    Predictions protein-protein interactions. Mostly wraps utils/ppi_prediction.py.
    """
    with open(output_ppis_fname, 'w+') as ppif, open(output_dpis_fname, 'w+') as dpif:
        writer = csv.writer(ppif, delimiter=',')
        
        for ui, uj in metadata.pairs:
            ppi_p, dpi_ps = utils.ppi_prediction.predict_ppi(ui, uj, metadata.domain_metadata, metadata.peptide_metadata, models)
            writer.writerow([ui,uj,ppi_p])
            dpif.write(json.dumps(dpi_ps) + "\n")

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--ppi_pairs", type=str, default=None)
    parser.add_argument("--domain_metadata", type=str, default="../data/predict/domain_metadata.csv")
    parser.add_argument("--peptide_metadata", type=str, default="../data/predict/peptide_metadata.csv")
    parser.add_argument("--output_ppi_prediction", type=str, default="ppi_predictions.csv",
            help="Output file for PPI predictions. Output is a csv file with format <ID 1>,<ID 2>,<PPI Likelihood>.")
    parser.add_argument("--output_dpi_prediction", type=str, default="dpi_predictions.txt",
            help="Output file for domain-peptide predictions. Output is a text file with each row representing a " + 
            "JSON string with domain-peptide interactions.")

    model_specification_group = parser.add_argument_group(title="Trained models specification")
    model_specification_group.add_argument("-m", "--models", type=str, default="models/hsm_pretrained/", 
            help="Defines a directory containing ." )
    model_specification_group.add_argument("--model_format", type=str, default="models/hsm_pretrained/models_format.csv",
            help="Defines the parameters for each model. Default: models/hsm_pretrained_format.csv")
    model_specification_group.add_argument("--amino_acids_order", type=str, default="models/amino_acids.txt"
            help="Define a different group (and order) of amino-acids. Needed to add a 'chemistry' type not currently used " + 
            "like phospho-serine/threonine. Default: models/amino_acids.txt.")
    
    opts = parser.parse_args()
    
    models_specification = utils.load_models.load_models_from_dir(opts.models, opts.model_format, opts.amino_acids_order)
    metadata_tuple = load_metadata(opts.domain_metadata, opts.peptide_metadata_fpath, models_specification, ppi_pairs_fpath=opts.ppi_pairs)
    predict_interactions(metadata_tuple, models_specification, opts.output_ppi_prediction, opts.output_dpi_prediction)
