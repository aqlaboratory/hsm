import os
import csv, json, pickle
from collections import *
from itertools import *

import numpy as np

from . import dpi_prediction

def _load_model(model_fpath, domain_len, peptide_len, is_fixed, amino_acid_ordering):
    # TODO: load models
    arr_archive = np.load(model_fpath)

    dw = arr_archive['domain_weights'].reshape(-1)
    pw = arr_archive['peptide_weights'].reshape(-1)

    if 'basis_weights' in arr_archive:
        iw = np.matmul(arr_archive['interaction_weights'], arr_archive['basis_weights']).reshape(-1)
    else:
        iw = arr_archive['interaction_weights'].reshape(-1) 
    b = arr_archive['bias']
    
    weights_tuple = dpi_prediction.DomainPeptideWeights(dw,pw,iw,b)

    model_fn = lambda domain_sequence, peptide_sequence: dpi_prediction.domain_peptide_interaction_prediction(domain_sequence, peptide_sequence,
            weights=weights_tuple, domain_length=domain_len, peptide_length=peptide_len, is_fixed=is_fixed, amino_acid_ordering=amino_acid_ordering)
    
    return model_fn
    

ModelsSpecificationFormat = namedtuple("ModelsSpecificationFormat", ["models_format", "models"]) 
def load_models_from_dir(models_dirpath, models_format_fpath, amino_acids_ordering):
    models_format = defaultdict(dict)
    
    model_filenames = dict()
    for domain_type, peptide_type, domain_align_len, peptide_align_len, is_fixed, model_fname in csv.reader(open(models_format_fpath, 'r')):
        domain_align_len, peptide_align_len = int(domain_align_len), int(peptide_align_len) 
        is_fixed = bool(int(is_fixed))
        
        models_format[domain_type][peptide_type] = (domain_align_len, peptide_align_len, is_fixed)
        
        model_filenames[(domain_type, peptide_type)] = model_fname 
    
    aa_ordering = {aa.strip():idx for idx, aa in enumerate(open(amino_acids_ordering, 'r'))}

    models = defaultdict(dict)
    for (domain_type, peptide_type), model_fname in model_filenames.items(): 
        fpath = os.path.join(models_dirpath, model_fname)
        assert os.path.exists(fpath)
        
        models[domain_type][peptide_type] = _load_model(fpath, *models_format[domain_type][peptide_type], aa_ordering)
    
    return ModelsSpecificationFormat(models_format, models)
