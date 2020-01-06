import os
from collections import *
from itertools import *

import numpy as np
import tensorflow as tf

from . import constants

binding_data = "binding.npy"
domain_sequences = "dseq_mtx.npy"
peptide_sequences = "pseq_mtx.npy"
interaction_sequences = "iseqs_mtx.npy"

InteractionDataTuple = namedtuple("InteractionDataTuple", ["binds", "domain_seqs", "peptide_seqs", "interaction_seqs"])
def split_data(data_directory, validation_chunk, include_all=False, excluded_chunks=None, n_folds=constants.DEF_n_folds, seed=constants.DEF_fold_seed):
    files = [constants.binding_file, constants.domain_file,  constants.peptide_file, constants.interaction_file]
    data = [np.load(os.path.join(data_directory, f)) for f in files]
    
    data_sz = data[0].shape[0]
    
    np.random.seed(seed)
    randomized = np.random.permutation(data_sz)
    chunks = np.array_split(randomized, n_folds)
    
    vindices = chunks[validation_chunk]
    
    validation_data = InteractionDataTuple(*[mtx[vindices] for mtx in data])

    if include_all:
        _excluded = []
    else:
        _excluded = [validation_chunk, *(excluded_chunks if excluded_chunks is not None else [])]
    
    tindices = [i for cidx, chunk in enumerate(chunks) for i in chunk if cidx not in _excluded]
    train_data = InteractionDataTuple(*[mtx[tindices] for mtx in data])

    np.random.seed()
    
    return train_data, validation_data

def _data_chunker(data, nchunks):
    randomized = np.random.permutation(data.binds.shape[0])
    splits = np.array_split(randomized, nchunks)

    for spl in splits:
        binds_spl = data.binds[spl]
        seqs_spl = [s[spl] for s in data[1:]]

        yield seqs_spl, binds_spl

def training_iterator(interaction_data, chunk_size=None, nchunks=None):
    if chunk_size is None and nchunks is None: 
        raise ValueError("Need to pass a chunk size or num. of chunks")

    if nchunks is None:
        nchunks = interaction_data.binds.shape[0] // chunk_size + 1

    while True:
        yield _data_chunker(interaction_data, nchunks) 

def make_sparse(sequences, domain_length, peptide_length, n_amino_acids):
    def _sparsify(arr, ncols):
        nrows = arr.shape[0]

        ridxes, compressed_cidxes = np.where(arr >= 0)
        cidxes = arr[ridxes, compressed_cidxes]
        
        vals = np.ones(ridxes.size)
        
        idxes = np.vstack([ridxes, cidxes]).T
         
        shape = [nrows, ncols]
        
        return tf.SparseTensorValue(
                    indices = idxes, 
                    values = vals,
                    dense_shape = shape
                )
    
    col_sizes = [domain_length * n_amino_acids, peptide_length * n_amino_acids, domain_length * peptide_length * n_amino_acids * n_amino_acids]
    return [_sparsify(m, ncols) for m, ncols in zip(sequences, col_sizes)]

