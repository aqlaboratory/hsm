import os
import json
from collections import *
from itertools import *
from tqdm import tqdm as tqdm

import numpy as np
import tensorflow as tf
from sklearn.metrics import roc_auc_score

from . import utils
from . import constants

class HSMDomainsModel(object):
    """
    Multi-domain model implementing the HSM / D model.
    
    Most arguments are passed into the model using the args initialization option. See constants.py for the 
    function for the keywords to use to set different options values.

    The required arguments are:
        - validation_chunk: an index that is validated against.
        - model_format: a ModelSpecificationTuple, as defined in constants.py
        - lambda_params: a nested python dictionary; nested for domain type then peptide type.
        - amino_acid_ordering: a dict that indexes the vectorized data.
        - args: all additional arguments; see constants.py
    """

    def __init__(self, validation_chunk, models_format, lambda_params, amino_acid_ordering, args=dict()):
        ## model_format inputs are namedtuples:
        # ModelSpecificationTuple = namedtuple("ModelSpecificationTuple", ["domain_type", "peptide_type", "domain_length", "peptide_length", "directory"])
        self.model_specs = dict()
        for model_format in models_format:
            self.model_specs[(model_format.domain_type, model_format.peptide_type)] = model_format
        self.n_amino_acids = len(amino_acid_ordering)

        # Set model hyper-parameters. Defaults to constants in the case where 
        self.epochs = args.get(constants.KEY_epochs, constants.DEF_epochs)
        self.data_split_seed = args.get(constants.KEY_chunk_seed, constants.DEF_fold_seed)
        self.n_folds = args.get(constants.KEY_n_folds, constants.DEF_n_folds)
        self.chunk_size = args.get(constants.KEY_chunk_size, constants.DEF_chunk_size)
        self.learning_rate = args.get(constants.KEY_learning_rate, constants.DEF_learning_rate)
        self.init_std = args.get(constants.KEY_standard_dev, constants.DEF_init_std)
        
        self.lambda_params = dict()
        for dtype, _lambdas in lambda_params.items():
            for ptype, lambdap in _lambdas.items():
                self.lambda_params[(dtype, ptype)] = lambdap
               
        self.validation_chunk = validation_chunk
        self.validate_step = args.get(constants.KEY_validate_step, constants.DEF_validate_step)
        self.exclude_indices = args.get(constants.KEY_exclude_chunks, None)
        self.include_all_data = args.get(constants.KEY_include_all_chunks, False)
        
        # Basis set size specification.
        self.basis_size = args.get(constants.KEY_basis_size, constants.DEF_basis_size)
        self.input_basis = args.get(constants.KEY_input_basis, None)
        self.train_basis = args.get(constants.KEY_train_basis, False)

        # Initialize tensorflow graph.
        self.sess = tf.Session()
        self.optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate) 
        
        self._variables(self.init_std)
        
        self._model(args)

        init_op = tf.group(tf.local_variables_initializer(), tf.global_variables_initializer())
        self.sess.run(init_op)

    
    def _variables(self, initialization_std):
        """
        Intializes variable. Shouldn't be called outside of the class (called from __init__).
        """
        # Define basis weights a bit earlier.
        basis_weights_shape = (self.basis_size, self.n_amino_acids ** 2)
        if self.input_basis is not None:
            # Load basis weights from pre-loaded option
            basis = np.load(self.input_basis)
            self.basis_weights = tf.Variable(basis, name='basis_weights', dtype=tf.float64, trainable=self.train_basis)
        else:
            # Randomly initialize basis weights
            self.basis_weights = tf.Variable(np.random.normal(scale=initialization_std, size=basis_weights_shape),
                    name="basis_weights", dtype=tf.float64, trainable=self.train_basis)
        
        PlaceholdersTuple = namedtuple("PlaceholdersTuple", ["domains", "peptides", "interacts", "binds"]) 
        WeightsTuple = namedtuple("WeightsTuple", ["domains", "peptides", "interacts", "bias"]) 

        # Placeholders and variables for different domain models 
        self.placeholders, self.weights = dict(), dict()
        for (dtype, ptype), model_fmt in self.model_specs.items():
            self.placeholders[(dtype, ptype)] = PlaceholdersTuple(
                    domains = tf.sparse.placeholder(dtype=tf.float64, name="domain_data_mtx_{0}_{1}".format(dtype, ptype)),
                    peptides = tf.sparse.placeholder(dtype=tf.float64, name="peptide_data_mtx_{0}_{1}".format(dtype, ptype)),
                    interacts = tf.sparse.placeholder(dtype=tf.float64, name="interaction_data_mtx_{0}_{1}".format(dtype, ptype)),
                    binds = tf.placeholder(dtype=tf.float64, name="binding_data_mtx_{0}_{1}".format(dtype, ptype))
                    )
            
            dlen, plen = model_fmt.domain_length, model_fmt.peptide_length 
            
            domain_weights_shape = (dlen * self.n_amino_acids, 1)
            dw = tf.Variable(np.random.normal(scale=initialization_std, size=domain_weights_shape), 
                    name="domain_weights_{0}_{1}".format(dtype, ptype), dtype=tf.float64)
            
            peptide_weights_shape = (plen * self.n_amino_acids, 1)
            pw = tf.Variable(np.random.normal(scale=initialization_std, size=peptide_weights_shape), 
                    name="peptide_weights_{0}_{1}".format(dtype, ptype), dtype=tf.float64)
          
            interact_weights_shape = (dlen * plen, self.basis_size)
            iw = tf.Variable(np.random.normal(scale=initialization_std, size=interact_weights_shape), 
                    name="interact_weights_{0}_{1}".format(dtype, ptype), dtype=tf.float64)
 
            bias = tf.Variable(0, name="bias_{0}_{1}".format(dtype, ptype), dtype=tf.float64)
            
            self.weights[(dtype, ptype)] = WeightsTuple(
                        domains=dw,
                        peptides=pw,
                        interacts=iw,
                        bias=bias
                    )

    def _model(self, args):
        """
        Defines model graph. Shouldn't be called outside of class (called from __init__).
        """
        # Different values can be accessed with these dictionaries throughout the model.
        self.sparsity_penalties = dict()
        self.logits, self.cross_entropies = dict(), dict()
        self.predict_ops = dict()

        for (dtype, ptype), model_fmt in self.model_specs.items():
            # Load placeholders, weights, and lambda for a specific model
            #   PlaceholdersTuple = namedtuple("PlaceholdersTuple", ["domains", "peptides", "interacts", "binds"]) 
            #   WeightsTuple = namedtuple("WeightsTuple", ["domains", "peptides", "interacts", "bias"]) 


            placeholders = self.placeholders[(dtype, ptype)]
            weights = self.weights[(dtype, ptype)]
            lambdap = self.lambda_params[(dtype, ptype)]
            
            # Scoping it for each domain-peptide model. Helpful to use w/ tensorboard.
            with tf.name_scope("compute_{0}_{1}".format(dtype, ptype)):
                # Sparsity Penalty. Combined across domains below.
                with tf.name_scope("sparsity"):
                    pen = (lambdap * tf.reduce_sum(tf.abs(weights.domains)) + 
                           lambdap * tf.reduce_sum(tf.abs(weights.peptides)) + 
                           lambdap * tf.reduce_sum(tf.abs(weights.interacts))) 
                    self.sparsity_penalties[(dtype, ptype)] = pen
                
                # Logit computation.
                with tf.name_scope("logits"):
                    # Basic model, computes logits from model weights.
                    domain_contrib = tf.sparse.sparse_dense_matmul(placeholders.domains, weights.domains)
                    peptide_contrib = tf.sparse.sparse_dense_matmul(placeholders.peptides, weights.peptides)
                    
                    _interact_weights = tf.expand_dims(tf.reshape(tf.matmul(weights.interacts, self.basis_weights), [-1]), -1)
                    interact_contrib = tf.sparse.sparse_dense_matmul(placeholders.interacts, _interact_weights)
                    
                    _logits = domain_contrib + peptide_contrib + interact_contrib + weights.bias            
                    self.logits[(dtype, ptype)] = _logits 
                
                # Cross entropy per-domain/peptide. Combined across domains below.
                self.cross_entropies[(dtype, ptype)] = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=placeholders.binds, logits=_logits, name="cross_entropy"))
                
                # Predict ops for each domain.
                self.predict_ops[(dtype, ptype)] = tf.math.sigmoid(_logits, name="predict") 
        
        # Loss Function computation and optimization.
        self.weight_penalty = tf.reduce_sum([v for v in self.sparsity_penalties.values()])
        self.cross_entropy = tf.reduce_sum([v for v in self.cross_entropies.values()])
        self.loss_function = self.cross_entropy + self.weight_penalty
        
        self.train_op = self.optimizer.minimize(self.loss_function)

    def optimize_step(self, data_iter_dct, nsteps):
        """
        A single training epoch. Passes in an iterator over data and 

        args:
            - data_iter_dct: a set of data iterators feeding in data of the form (sequence_tuple, binds_arr)
                sequence_tuple is ordered[domain sequences, peptide sequences, interact sequences]                
        """
        costs = list() 
        
        for step in range(nsteps):
            feed_dict = dict()
                
            for (dtype, ptype), model_format in self.model_specs.items():
                # PlaceholdersTuple = namedtuple("PlaceholdersTuple", ["domains", "peptides", "interacts", "binds"])
                placeholders = self.placeholders[(dtype, ptype)]
                
                # Convert to sparse tensors
                sequences, binds = next(data_iter_dct[(dtype, ptype)])
                sparse_seqs = utils.make_sparse(sequences, model_format.domain_length, model_format.peptide_length, self.n_amino_acids)
                
                # _p is the placeholder matched with the value, _v
                for _p, _v in zip(placeholders, [*sparse_seqs, np.expand_dims(binds, -1)]):
                    feed_dict[_p] = _v
            
            train_cost, _ = self.sess.run([self.loss_function, self.train_op], feed_dict=feed_dict)
            costs.append(train_cost)
        
        return costs
           
    def predict(self, data_dct):
        """
        Make a single prediction. 

        args:
            - data_dct: a data tuple of the form InteractionDataTuple as defined in utils.py 
        """
        predictions_dct = dict()
        
        # This can be done in one call. Splitting it up for reader ease.
        for (dtype, ptype), data in data_dct.items():
            binds, sequences = data[0], list(data[1:])
        
            # Convert to sparse tensors.
            model_fmt = self.model_specs[(dtype, ptype)]
            dseqs_sp, pseqs_sp, iseqs_sp = utils.make_sparse(sequences, model_fmt.domain_length, model_fmt.peptide_length, self.n_amino_acids)
            
            placeholder = self.placeholders[(dtype, ptype)]

            feed_dict = {
                    placeholder.domains: dseqs_sp,
                    placeholder.peptides: pseqs_sp,
                    placeholder.interacts: iseqs_sp,
                    placeholder.binds: np.expand_dims(binds, -1)
                    }
            
            # Get domain-specific predict_op and predict.
            predict_op = self.predict_ops[(dtype, ptype)]
            predictions = self.sess.run(predict_op, feed_dict=feed_dict)

            """
            # Realized this is unneccessary debugging -- was running on
            # an unstable GPU (2 on exx1)
            print(predictions)
            nanpreds = np.isnan(predictions)
            if sum(nanpreds) >= 1:
               print("Predictions contain nan", predictions[nanpreds])
            posinfpreds = np.isposinf(predictions)
            if sum(posinfpred) >= 1:
               print("Predictions contain +inf", predictions[posinfpreds])
            neginfpreds = np.isneginf(predictions)
            if sum(neginfpreds) >= 1:
               print("Predictions contain -inf", predictions[neginfpreds])
            """

            roc_auc = roc_auc_score(binds, predictions) 
            
            predictions_dct[(dtype, ptype)] = (binds, predictions, roc_auc)
        
        return predictions_dct 
    
    def train(self, data_directory):
        """
        Fully train the model given the input parameters. Runs for the input / default number of epochs. 
        Mostly just wraps the optimize_step function. 

        args:
            - data_directory: the input data directory. The data is split using the input metadata.  
        """

        # Initialize train and test data splits.
        train_data_dct, test_data_dct = dict(), dict() 
        for (dtype, ptype), model_fmt in self.model_specs.items():
            # model_fmt is a ModelSpecificationTuple = namedtuple("ModelSpecificationTuple", ["domain_type", "peptide_type", "domain_length", "peptide_length", "directory"])
            domain_data_directory = os.path.join(data_directory, model_fmt.directory)
            train_data, test_data = utils.split_data(domain_data_directory, self.validation_chunk, include_all=self.include_all_data, excluded_chunks=self.exclude_indices, n_folds=self.n_folds, seed=self.data_split_seed)
            
            train_data_dct[(dtype, ptype)] = train_data
            test_data_dct[(dtype, ptype)] = test_data
        
        # Variable input sizes. Calculate the maximum number of chunks needed for any one domain
        #   given the input max chunk size.
        n_chunks = max((v[0].shape[0] // self.chunk_size) + 1 for v in train_data_dct.values())
        
        # Initial data iterators.
        data_iterators_dct = dict()
        for (dtype, ptype), train_data in train_data_dct.items():
            data_iterators_dct[(dtype, ptype)] = utils.training_iterator(train_data, nchunks=n_chunks)

        self.costs = list() 
        self.aucs = list()
        self.aucs_per_dp = list()
        for epoch in tqdm(range(self.epochs), desc="Epochs"):
            epoch_iters_dct = {k:next(v) for k,v in data_iterators_dct.items()}
            epoch_costs = self.optimize_step(epoch_iters_dct, n_chunks)
            self.costs.append(epoch_costs)
            
            if (epoch + 1) % self.validate_step == 0:
                _perf = self.predict(test_data_dct)
                perf_dct = dict()
                for (dtype, ptype), (_, _, auc) in _perf.items():
                    print("Epoch: {0} (AUC: {1}; Domain: {2}, Peptide: {3})".format(epoch+1, auc, dtype, ptype))
                    perf_dct['%s_%s' % (dtype, ptype)] = auc
                    self.aucs.append(auc)
                self.aucs_per_dp.append(perf_dct)

        self.final_predictions = self.predict(test_data_dct)

    def save_model(self, output_directory, label=None):
        """
        Save model outputs from a trained model. The output files are of the form: <id>.<optional_label>.<spec>.<ext>
        The id is defined as the current time (as a string). An optional label istring can be supplied.
        The spec defines the output, either metadata or the domain type of output. ext is the file
        extension, either a numpy file or a json file.

        args:
            - output_directory: a string defining a directory to output the saved model to.  
        """
        import time        
        
        model_output_id = str(time.time()).replace("-", "_")
        if label is not None: model_output_id = model_output_id+'.'+label
        
        # Metadata. Saved as dict (to JSON file) with four keys defining model type, 
        # model format, the results of the last prediction, and parameters
        results = list()
        for (dtype, ptype), preds in self.final_predictions.items():
            _r = {"binds": list(map(int, preds[0])),
                  "predictions": list(map(float, preds[1])),
                  "auc": float(preds[2])
                  }

            results.append([[dtype, ptype], _r]) 
        output = {
                "Model Type": "HSM / ID",
                "Model Formats": [list(fmt) for fmt in self.model_specs.values()], 
                "results": results, 
                "parameters": {
                        "costs": self.costs,
                        "validation chunk": self.validation_chunk,
                        constants.KEY_chunk_seed: self.data_split_seed,
                        constants.KEY_n_folds: self.n_folds,
                        constants.KEY_epochs: self.epochs,
                        constants.KEY_learning_rate: self.learning_rate
                    }
                }
        
        ometadata_file = os.path.join(output_directory, "{0}.metadata.json".format(model_output_id))
        with open(ometadata_file, 'w+') as f:
            json.dump(output, f)
        
        # The model, saved as an array archive (.npz) file.
        for (dtype, ptype), weights in self.weights.items():
            model_file = os.path.join(output_directory, "{0}.{1}.{2}.npz".format(model_output_id, dtype, ptype))
            
            # WeightsTuple = namedtuple("WeightsTuple", ["domains", "peptides", "interacts", "bias"])
            arr = [weights.domains, weights.peptides, weights.interacts, self.basis_weights, weights.bias]
            dw, pw, iw, basis, bias = self.sess.run(arr)
            np.savez(model_file, domain_weights=dw, peptide_weights=pw, interact_weights=iw, basis_weights=basis, bias=bias)
