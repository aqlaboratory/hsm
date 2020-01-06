"""
This file provides a CLI for running models with a variety of possible 
options.

Note: you can also import the models directly and run them.

The easiest way to understand these config options is likely by running the
built-in help method, i.e.:
    $ python <path_to_dir>/train.py --help
"""
import os
import csv
import argparse

from collections import *
from itertools import *

import models
from models import constants as constants

def create_parser():
    parser = argparse.ArgumentParser(description="A wrapper to specify, train, and assess different models. " +
            "Models are implemented primarily in the models/ directory and may be accessed as a package."
            )
    
    parser.add_argument("validation_chunk_idx", type=int, default=0, help="Validation chunk")

    model_spec_group = parser.add_argument_group(title="Model Specification")
    
    # Specificy model type 
    model_group = model_spec_group.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--hsm_id", action="store_true", default=False, 
            help="Train the HSM/ID model. Domain specification must specify a single domain." ) 
    model_group.add_argument("--hsm_d", action="store_true",  
            help="Train the HSM/D model. Domain specification may include one or multiple domains.")
    model_spec_group.add_argument("-d", "--domains", nargs='+', type=str,
            help="Specify a domain type for training a new model. To specify a specific domain - peptide type, pass with a ':' character separating types.")
    model_spec_group.add_argument("-a", "--include_all_domains", action='store_true', default=False,
            help="Include all possible domains for which data are available")
    
    hyperparam_group = parser.add_argument_group(title="Learning hyperparameters")
    hyperparam_group.add_argument("-l", "--lambdap", type=float, default=constants.DEF_lambda, 
            help="Sparsity parametrization parameter. Default: {0}".format(str(constants.DEF_lambda)))
    hyperparam_group.add_argument("--lambda_params", type=str, nargs='+',
            help="Sparsity parametrization per-domain, passed as a string of format {domain}:{lambda} or {domain}:{peptide}:{lambda} where each choice" + 
            "is comma-separated. Default: {0}".format(str(constants.DEF_lambda)))
    hyperparam_group.add_argument("-r", "--learning_rate", type=float, default=constants.DEF_learning_rate, 
            help="Add an initial learning rate. Default: {0}".format(str(constants.DEF_learning_rate)))
    hyperparam_group.add_argument("-s", "--standard_deviation", type=float, default=constants.DEF_init_std, 
            help="Set the initialization standard dev. Default: {0}".format(str(constants.DEF_init_std)))
    hyperparam_group.add_argument("-e", "--epochs", type=int, default=constants.DEF_epochs, 
            help="Max number of epochs to run the strategy for. Default: {0}".format(str(constants.DEF_epochs)))
    hyperparam_group.add_argument("-v", "--validate_step", type=int, default=constants.DEF_validate_step,
            help="The step number to run validation on. Default: {0}.".format(str(constants.DEF_validate_step)))
    
    data_handling_group = parser.add_argument_group(title="Data handling group") 
    data_handling_group.add_argument("--chunk_seed", type=int, default=0,
            help="Seed used to synchronize chunks across a single split. Use the same seed to train on different chunks.")
    data_handling_group.add_argument("--exclude_chunks", nargs="+", type=int, 
            help="Exclude specific chunks from training or testing. Helpful for excluding a test chunk")
    data_handling_group.add_argument("--include_all_data", action='store_true', default=False, 
            help='Use all data in training. Used at end to get final predictions.')
    data_handling_group.add_argument("-c", "--chunk_size", type=int, default=constants.DEF_chunk_size,
            help="The maximum chunk size for any iteration. Actual chunk size is less than this (calculated as the " +
            "amount of data divided by chunk size plus one. Default: {0}.".format(str(constants.DEF_chunk_size)))
    data_handling_group.add_argument("--n_folds", type=int, default=constants.DEF_n_folds, 
            help="Number of chunks to split dataset into. Default: {0}".format(constants.DEF_n_folds))

    # Input / output data processing.
    input_output_group = parser.add_argument_group(title="Input/output group")
    input_output_group.add_argument("-i", "--input_directory", type=str, default="../data/train/hsm_preprocessed/",
            help="Input data directory. The form of this directory should be the same as output by convert_binding_data.py")
    input_output_group.add_argument("-o", "--output_directory", type=str, default="outputs/", 
            help="Location to output data to.")
    input_output_group.add_argument("--amino_acids_ordering", type=str, default="../data/amino_acid_ordering.txt",
            help="Amino acid ordering file. Allows user to add new amino acids (e.g. phospho-ser/thr).")

    # Model specific arguments:
    basis_spec_group = parser.add_argument_group(title="Basis Specification")
    basis_spec_group.add_argument("-b", "--basis_size", type=int, default=constants.DEF_basis_size,
            help="Set the hierarchical/shared-basis size. Default: 100".format(str(constants.DEF_basis_size)))
    basis_spec_group.add_argument("--load_basis", type=str, 
            help="Run the model w/ a pre-defined basis. Basis is passed as a numpy file.")
    basis_spec_group.add_argument("--no_train", action="store_true", default=False,
            help="Makes the basis set non-trainable for a given initialization.")
    
    return parser

def format_arguments_directory(options):
    arguments = dict()
    
    arguments[constants.KEY_learning_rate] = options.learning_rate
    arguments[constants.KEY_standard_dev] = options.standard_deviation
    arguments[constants.KEY_epochs] = options.epochs
    arguments[constants.KEY_validate_step] = options.validate_step
    
    arguments[constants.KEY_chunk_seed] = options.chunk_seed
    arguments[constants.KEY_exclude_chunks] = options.exclude_chunks
    arguments[constants.KEY_include_all_chunks] = options.include_all_data
    arguments[constants.KEY_chunk_size] = options.chunk_size

    arguments[constants.KEY_basis_size] = options.basis_size
    arguments[constants.KEY_input_basis] = options.load_basis
    arguments[constants.KEY_train_basis] = options.no_train
    
    arguments[constants.KEY_output_directory] = options.output_directory
    arguments[constants.KEY_n_folds] = options.n_folds

    
    arguments = {k:v for k,v in arguments.items() if v is not None}
    
    return arguments
    
    
def get_model(args, opts):
    def _match_model(dinput, model_spec_dir):
        if ":" in dinput:
            dtype, ptype = dinput.split(":")
        else:
            dtype = dinput
            ptype = None

        assert dtype in possible_models 
        
        if ptype is None:
            assert len(possible_models[dtype]) == 1
            return next(iter(possible_models[dtype].values()))
        else:
            assert ptype in possible_models[dtype]
            return possible_models[dtype][ptype]

    model_formats_ifile = os.path.join(opts.input_directory, constants.model_specification_file)
    
    all_models = list()
    possible_models = defaultdict(dict)
    for dtype, ptype, dlen, plen, idir in csv.reader(open(model_formats_ifile, 'r'), delimiter=','):
        mformat = constants.ModelSpecificationTuple(dtype, ptype, int(dlen), int(plen), idir)
        possible_models[dtype][ptype] = mformat
        all_models.append(mformat)
    
    def _process_lambdas(lambda_param, lambda_params, possible_models=possible_models):
        lambda_dirs = defaultdict(dict)
        if lambda_params is None:
            for dtype, model_types in possible_models.items():
                lambda_dirs[dtype] = {ptype:lambda_param for ptype in possible_models[dtype].keys()}
        else:
            for arg in lambda_params:
                spl = arg.split(":")
                dtype, lp = spl[0], float(spl[-1])

                if len(spl) == 2:
                    lambda_dirs[dtype] = {ptype:lp for ptype in possible_models[dtype].keys()}
                elif len(spl) == 3:
                    ptype = spl[1]
                    lambda_dirs[dtype][ptype] = lp
                else:
                    raise ValueError("Lambdas mis-specified")
            
            for dtype, model_types in possible_models.items():
                for ptype in model_types.keys():
                    if ptype not in lambda_dirs[dtype]:
                        lambda_dirs[dtype][ptype] = constants.DEF_lambda
        
        return lambda_dirs
    
    lambda_params_dct = _process_lambdas(options.lambdap, options.lambda_params)

    if opts.hsm_id:
        assert len(opts.domains) == 1
        
        model_formats = _match_model(opts.domains[0], possible_models)
        model_type = models.hsm_id.HSMIndependentDomainsModel
        input_directories = os.path.join(opts.input_directory, model_formats.directory) 
    elif len(opts.domains) == 1:
        model_formats = _match_model(opts.domains[0], possible_models)
        model_type = models.hsm_d_singledomain.HSMSingleDomainsModel
        input_directories = os.path.join(opts.input_directory, model_formats.directory) 
    elif opts.include_all_domains:
        model_formats =  model_formats_ifile
        model_type = models.hsm_d.HSMDomainsModel
        input_directories = opts.input_directory
    else:
        model_formats = [_match_model(d, possible_models) for d in opts.domains]
        model_type = models.hsm_d.HSMDomainsModel
        input_directories = opts.input_directory
    
    amino_acid_ordering = {aa.strip():idx for idx, aa in enumerate(open(opts.amino_acids_ordering, 'r'))}

    model = model_type(opts.validation_chunk_idx, model_formats, lambda_params_dct, amino_acid_ordering, args) 
    return model, input_directories

if __name__=='__main__':
    parser = create_parser()
    options = parser.parse_args()
    assert os.path.exists(options.input_directory) and os.path.exists(options.output_directory) and os.path.exists(options.amino_acids_ordering)

    arguments = format_arguments_directory(options)

    model, idirs = get_model(arguments, options)
    model.train(idirs)
    model.save_model(options.output_directory)
