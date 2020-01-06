import os
import csv
from collections import *
from itertools import *
from tqdm import tqdm as tqdm

import numpy as np

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
    
    # Pad to equal length
    npad = len(sequence)  - len(vectorized)
    vectorized.extend(-1 for _ in range(npad))

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
    
    # Pad to equal length
    npad = len(domain_sequence) * len(peptide_sequence) - len(vectorized)
    vectorized.extend(-1 for _ in range(npad))

    return vectorized

def _load_binding_data(input_data_file, amino_acid_ordering,  progressbar=False):
    """
    Input data format should be a csv file of the form:
        Domain-Type,Domain-Protein-Identifier,Aligned-Domain-Sequence,Peptide-Type,Peptide-Protein-Identifier,Aligned-Peptidic-Sequence

    An iterator that returns data grouped by model type. Model types are automatically inferred from the input data file.
    """
    get_model_type = lambda row: (row[0], row[2], len(row[1]), len(row[3]))
    
    model_types, total = set(), 0
    for row in csv.reader(open(input_data_file, 'r')):
        total += 1
        model_type = get_model_type(row)

        model_types.add(model_type)

    for model_type in tqdm(model_types, disable=(not progressbar), desc="Model types"):
        binds = list()
        domain_seqs, peptide_seqs, interact_seqs = list(), list(), list()

        for row in tqdm(csv.reader(open(input_data_file, 'r')), disable=(not progressbar), desc="Data processing", total=total):
            if get_model_type(row) != model_type: continue

            b = 1 if float(row[4]) > 0 else 0 
            binds.append(b)

            domain_seqs.append(_vectorize_sequence(row[1], amino_acid_ordering))
            peptide_seqs.append(_vectorize_sequence(row[3], amino_acid_ordering))
            interact_seqs.append(_vectorize_interaction(row[1], row[3], amino_acid_ordering))
        
        binds = np.array(binds)
        vectorized = [np.array(a) for a in [domain_seqs, peptide_seqs, interact_seqs]]

        yield model_type, vectorized, binds

def convert_binding_data(input_data_file, output_data_directory, amino_acid_ordering, progressbar=False):
    """
    Function that converts data. Mostly wraps functions above.
    """

    assert os.path.exists(output_data_directory)
    
    model_fmt = list()

    processed_binding_data = defaultdict(list)
    for model_type, seqs_vectorized, binds in _load_binding_data(input_data_file, amino_acid_ordering, progressbar=progressbar):
        model_odirname = "{0}_{1}".format(*model_type)
        model_odirpath = os.path.join(output_data_directory,model_odirname) 
        os.mkdir(model_odirpath)
        
        model_fmt.append((*model_type, model_odirname))
        
        np.save(os.path.join(model_odirpath, 'binding.npy'), np.array(list(binds)))
        
        output_files = ["dseq_mtx.npy", "pseq_mtx.npy", "iseqs_mtx.npy"]
        for seqs, ofname in zip(seqs_vectorized, output_files):
            np.save(os.path.join(model_odirpath, ofname), seqs)
        
    with open(os.path.join(output_data_directory, 'amino_acid_ordering.txt'), 'w+') as f:
        amino_acids_list = [aa for aa, idx in sorted(amino_acid_ordering.items(), key=lambda x: x[1])]
        f.write('\n'.join(amino_acids_list))
    
    with open(os.path.join(output_data_directory, 'models_specification.csv'), 'w+') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(model_fmt)
   
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_data_file", type=str)
    parser.add_argument("output_data_directory", type=str)
    parser.add_argument("-a", "--amino_acid_ordering", type=str, default="../data/amino_acid_ordering.txt")
    parser.add_argument("-p", "--progressbar", action='store_true', default=False)
    opts = parser.parse_args()

    aa_order = {aa.strip():idx for idx, aa in enumerate(open(opts.amino_acid_ordering, 'r'))}
    convert_binding_data(opts.input_data_file, opts.output_data_directory, aa_order, progressbar=opts.progressbar)
