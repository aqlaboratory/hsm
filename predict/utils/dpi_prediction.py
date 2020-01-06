from collections import *
from itertools import *
import numpy as np

DomainPeptideWeights = namedtuple("DomainPeptideWeights", ["domain_weights", "peptide_weights", "interaction_weights", "bias"])

def _vectorize_sequence(sequence, amino_acid_ordering):
     """
     Computes a one-hot embedding of an input sequence.  
     
     Returns:
         - list. Non-zero indices of one-hot embedding matrix of a sequence.
                 Non-flattened, this matrix has dimensions:
                     (sequence length, # of amino acids represented)
     """
     aa_len = len(amino_acid_ordering)
 
     vectorized = list()
 
     sequence_indexed = [(sidx, saa) for sidx, saa in enumerate(sequence) if saa in amino_acid_ordering]
     for sidx, saa in sequence_indexed:
         idxed = sidx * aa_len + amino_acid_ordering[saa]
 
         vectorized.append(idxed)
 
     return vectorized

def _vectorize_interaction(domain_sequence, peptide_sequence, amino_acid_ordering):
    """
    Computes a one-hot embedding of the interaction between the domain- and peptidic-sequence.
    
    Returns: 
        - list. Non-zero indices for the interaction between domain and peptide sequences. 
                Non-flattened, this matrix has dimensions:
                    (domain sequence length, peptide sequence length, # of amino acids represented, # of amino acids represented)

    """

    aa_len = len(amino_acid_ordering)
    domain_idx_offset = len(peptide_sequence) * aa_len * aa_len
    peptide_idx_offset = aa_len * aa_len

    vectorized = list()

    domain_indexed = [(didx, daa) for didx, daa in enumerate(domain_sequence) if daa in amino_acid_ordering]
    peptide_indexed = [(pidx, paa) for pidx, paa in enumerate(peptide_sequence) if paa in amino_acid_ordering]

    for (didx, daa), (pidx, paa) in product(domain_indexed, peptide_indexed):

        idxed = didx * domain_idx_offset + pidx * peptide_idx_offset + amino_acid_ordering[daa] * aa_len + amino_acid_ordering[paa]
        vectorized.append(idxed)

    return vectorized

def _compute_interaction(dseq, pseq, weights, amino_acid_ordering):
    """
    Compute the likelihood of a given domain and peptide sequence given an input set of weights.
    Assumes domain and peptide sequence are aligned.

    args:
        - dseq: domain sequence as string
        - pseq: peptide sequence as string
        - weights: input interaction weights
        - amino_acid_ordering: ordering of amino acids
    """

    domain_contrib = sum([weights.domain_weights[didx] for didx in _vectorize_sequence(dseq, amino_acid_ordering)])
    peptide_contrib = sum([weights.peptide_weights[pidx] for pidx in _vectorize_sequence(pseq, amino_acid_ordering)])
    interact_contrib = sum([weights.interaction_weights[iidx] for iidx in _vectorize_interaction(dseq, pseq, amino_acid_ordering)])
    
    logit = domain_contrib + peptide_contrib + interact_contrib + weights.bias

    return 1 / (1 + np.exp(-logit))

def emit_peptide_sequences(pseq, plen, overhang=4):
    """
    Scans over variable length peptide region and outputs sequences.
    """
    if len(pseq) < plen:
        diff = plen - len(pseq)
        for roffs in range(0,diff+1):
            loffs = diff - roffs
            yield loffs * '-' + pseq + roffs * '-'
    elif len(pseq) > plen:
        for i in range(0, (len(pseq) + 1) - plen):
            yield pseq[i:i+plen]
    else:
        yield pseq

def _slide_interactions(dseq, pseq, weights, peptide_length, amino_acid_ordering):
    """
    Computes the likelihood of a sliding-window interaction.
    """

    likelihood_dist = np.array([_compute_interaction(dseq, pseq_subseq, weights, amino_acid_ordering) for pseq_subseq in emit_peptide_sequences(pseq, peptide_lengtht)])
    
    no_bind = np.prod(1-likelihood_dist)

    bind = sum(no_bind * (i / (1-i)) for i in likelihood_dist)
    return bind / (bind + no_bind)

def domain_peptide_interaction_prediction(domain_sequence, peptide_sequence, weights, domain_length, peptide_length, is_fixed, amino_acid_ordering):
    if is_fixed:
        return _compute_interaction(domain_sequence, pseq, weights, amino_acid_ordering)
    else:
        return _slide_interactions(domain_sequence, peptide_sequence, weights, peptide_length, amino_acid_ordering)
