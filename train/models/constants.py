

"""
Hyper-parameters for model training
"""
DEF_lambda = 1e-5
DEF_learning_rate = 1e-4
DEF_init_std = 1e-2
DEF_epochs = 100
DEF_validate_step = 100
DEF_chunk_size = 512

DEF_basis_size = 100

DEF_fold_seed = 0
DEF_n_folds = 8

KEY_validation_chunk = "Validation Chunk Index"

KEY_learning_rate = "Learning Rate"
KEY_lambdas = "Lambda Params"
KEY_standard_dev = "Standard Deviation"
KEY_epochs = "Epochs"
KEY_validate_step = "Validate Step"

KEY_chunk_seed = "Chunking Seed" 
KEY_exclude_chunks = "Exclude Chunks"
KEY_include_all_chunks = "Include All Chunks"
KEY_chunk_size = "Chunk Size"
KEY_n_folds = "Number of folds"

KEY_basis_size = "Basis Size"
KEY_input_basis = "Input Basis Filepath"
KEY_train_basis = "Train Basis"

KEY_output_directory = "Output Directory"

binding_file = "binding.npy"
domain_file = "dseq_mtx.npy"
peptide_file = "pseq_mtx.npy"
interaction_file = "iseqs_mtx.npy"

model_specification_file = "models_specification.csv"

from collections import namedtuple
ModelSpecificationTuple = namedtuple("ModelSpecificationTuple", [
    "domain_type", 
    "peptide_type", 
    "domain_length", 
    "peptide_length", 
    "directory"
    ])
