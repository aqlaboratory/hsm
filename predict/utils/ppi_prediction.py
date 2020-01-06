from collections import *
from itertools import *
import numpy as np

from .timeout_error import timeout, TimeoutError

def is_valid(i, j, domain_metadata, peptide_metadata, models):
    di, pi = domain_metadata.get(i, list()), peptide_metadata.get(i, list())
    dj, pj = domain_metadata.get(j, list()), peptide_metadata.get(j, list())
    
    def _is_valid_oneway(domains, peptides):
        dtypes = set(d[0] for d in domains)
        valid_ptypes = set(p for d in dtypes for p in models.models[d].keys())

        ptypes = set(p[0] for p in peptides)
        return len(ptypes.intersection(valid_ptypes)) > 0
    
    return _is_valid_oneway(di,pj) or _is_valid_oneway(dj,pi)

@timeout(2)
def _predict_ppi(interactions, max_ensemble_size):
    def _valid_ensembles(ensemble_size, interactions=interactions):
        for idxes in combinations(range(len(interactions)), r=ensemble_size):
            peptides = set(tuple(interactions[i][1]) for i in idxes)
            domains = set(tuple(interactions[i][0]) for i in idxes)
            
            if len(peptides) == len(domains) == ensemble_size:
                yield list(idxes)
    
    binding = np.array([i[-1] for i in interactions])

    not_bound = np.prod(1-binding)

    bound = 0
    for ensemble_size in range(1, max_ensemble_size + 1):
        changed = False

        for ensemble_idxes in _valid_ensembles(ensemble_size):
            factor = np.prod(binding[ensemble_idxes]) / (np.prod(1 - binding[ensemble_idxes]))
            bound += not_bound * factor
            changed = True

        if not changed: break
    
    return bound / (not_bound + bound)

def predict_ppi(i, j, domain_metadata, peptide_metadata, models):
    potential_interactions = chain(product(domain_metadata.get(i,list()), peptide_metadata.get(j, list())),
                                   product(domain_metadata.get(j,list()), peptide_metadata.get(i, list())))
    
    interactions = list()
    for (dtype, dseq), (ptype, pseq) in potential_interactions:
        if dtype in models.models and ptype in models.models[dtype]:
            dpi_p = models.models[dtype][ptype](dseq, pseq)
            interactions.append([(dtype, dseq), (ptype, pseq), dpi_p])

    if len(interactions) == 0:
        return "No Interactions", interactions
    
    max_ensemble_size = min(len(set(i[0] for i in interactions)), len(set(i[1] for i in interactions)))
    
    try:
        ppi_p = _predict_ppi(interactions, max_ensemble_size)
        return ppi_p, [[list(dtup), list(ptup), v] for dtup, ptup, v in interactions]
    except TimeoutError:
        return "Timeout Error, increase compute time", interactions
