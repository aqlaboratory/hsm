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

class HSMIndependentDomainsModel(object):
    """
    Single-domain model implementing the HSM / ID model.
    
    Most arguments are passed into the model using the args initialization option. See constants.py for the 
    function for the keywords to use to set different options values.

    The required arguments are:
        - validation_chunk: an index that is validated against.
        - model_format: a ModelSpecificationTuple, as defined in constants.py
        - lambda_params: a nested python dictionary; nested for domain type then peptide type.
        - amino_acid_ordering: a dict that indexes the vectorized data.
        - args: all additional arguments; see constants.py
    """

    def __init__(self, validation_chunk, model_format, lambda_params, amino_acid_ordering, args=dict()):
        ## model_format inputs are namedtuples:
        # ModelSpecificationTuple = namedtuple("ModelSpecificationTuple", ["domain_type", "peptide_type", "domain_length", "peptide_length", "directory"])
        self.model_spec = model_format
        self.domain_length = model_format.domain_length
        self.peptide_length = model_format.peptide_length
        self.n_amino_acids = len(amino_acid_ordering)
        
        # Set model hyper-parameters. Defaults to constants in the case where 
        self.epochs = args.get(constants.KEY_epochs, constants.DEF_epochs)
        self.data_split_seed = args.get(constants.KEY_chunk_seed, constants.DEF_fold_seed)
        self.n_folds = args.get(constants.KEY_n_folds, constants.DEF_n_folds)
        self.chunk_size = args.get(constants.KEY_chunk_size, constants.DEF_chunk_size)
        self.learning_rate = args.get(constants.KEY_learning_rate, constants.DEF_learning_rate)
        self.init_std = args.get(constants.KEY_standard_dev, constants.DEF_init_std)
        self.lambda_param = lambda_params[model_format.domain_type][model_format.peptide_type] 
        
        self.validation_chunk = validation_chunk
        self.validate_step = args.get(constants.KEY_validate_step, constants.DEF_validate_step)
        self.exclude_indices = args.get(constants.KEY_exclude_chunks, None)
        self.include_all_data = args.get(constants.KEY_include_all_chunks, False)

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
        # Placeholders for different input data.
        self.domain_data = tf.sparse.placeholder(dtype=tf.float64, name="domain_data_mtx")
        self.peptide_data = tf.sparse.placeholder(dtype=tf.float64, name="peptide_data_mtx")
        self.interact_data = tf.sparse.placeholder(dtype=tf.float64, name="interaction_data_mtx")
        self.binding_data = tf.placeholder(dtype=tf.float64, name="binding_data_mtx")
        
        # Weights variables to be trained in the model.
        domain_weights_shape = (self.domain_length * self.n_amino_acids, 1)
        self.domain_weights = tf.Variable(np.random.normal(scale=initialization_std, size=domain_weights_shape),
                name="domain_weights", dtype=tf.float64)

        peptide_weights_shape = (self.peptide_length * self.n_amino_acids, 1)
        self.peptide_weights = tf.Variable(np.random.normal(scale=initialization_std, size=peptide_weights_shape),
                name="peptide_weights", dtype=tf.float64)

        interaction_weights_shape = (self.domain_length * self.peptide_length * (self.n_amino_acids ** 2), 1)
        self.interaction_weights = tf.Variable(np.random.normal(scale=initialization_std, size=interaction_weights_shape), 
                name="interaction_weights", dtype=tf.float64)

        self.bias = tf.Variable(0, name="bias", dtype=tf.float64)

    def _model(self, args):
        """
        Defines model graph. Shouldn't be called outside of class (called from __init__).
        """
        with tf.name_scope("sparsity"):
            # L1-sparsity penalty.
            self.weight_penalty = (self.lambda_param * tf.reduce_sum(tf.abs(self.domain_weights)) + 
                                  self.lambda_param * tf.reduce_sum(tf.abs(self.peptide_weights)) + 
                                  self.lambda_param * tf.reduce_sum(tf.abs(self.interaction_weights))) 
            
        with tf.name_scope("logits"):
            # Basic model, computes logits from model weights.
            domain_contrib = tf.sparse.sparse_dense_matmul(self.domain_data, self.domain_weights) 
            peptide_contrib = tf.sparse.sparse_dense_matmul(self.peptide_data, self.peptide_weights)
            interact_contrib = tf.sparse.sparse_dense_matmul(self.interact_data, self.interaction_weights)

            self.logits = domain_contrib + peptide_contrib + interact_contrib + self.bias
        
        # Loss Function computation and optimization.
        self.cross_entropy = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=self.binding_data, logits=self.logits, name="cross_entropy"))
        self.loss_function = self.cross_entropy + self.weight_penalty
        
        self.train_op = self.optimizer.minimize(self.loss_function)
        
        # Ops for computing predictions.
        self.predict_op = tf.math.sigmoid(self.logits, name="predict_sigmoid")

    def optimize_step(self, data_iter):
        """
        A single training epoch. Passes in an iterator over data and 

        args:
            - data_iter: a data iterator feeding in data of the form (sequence_tuple, binds_arr)
                sequence_tuple is ordered[domain sequences, peptide sequences, interact sequences]                
        """
        costs = list() 
        for sequences, binds in data_iter:
            # Convert to sparse tensors for use in the model.
            dseqs_sp, pseqs_sp, iseqs_sp = utils.make_sparse(sequences, self.domain_length, self.peptide_length, self.n_amino_acids)
            
            # Using feed_dicts is a bit unfortunate. It cleans up the data file handling on the
            # front end. To optimize, switch over to TensorFlow's data processing pipeline.
            feed_dict = {
                    self.domain_data: dseqs_sp,
                    self.peptide_data: pseqs_sp,
                    self.interact_data: iseqs_sp,
                    self.binding_data: np.expand_dims(binds, -1)
                    }
            
            # Compute step.
            train_cost, _ = self.sess.run([self.loss_function, self.train_op], feed_dict=feed_dict)
            costs.append(train_cost)

        return costs
    
    def predict(self, data):
        """
        Make a single prediction. 

        args:
            - data: a data tuple of the form InteractionDataTuple as defined in utils.py 
        """
        binds, sequences = data[0], list(data[1:])
        
        # Convert to sparse tensors.
        dseqs_sp, pseqs_sp, iseqs_sp = utils.make_sparse(sequences, self.domain_length, self.peptide_length, self.n_amino_acids)

        feed_dict = {
                self.domain_data: dseqs_sp,
                self.peptide_data: pseqs_sp,
                self.interact_data: iseqs_sp,
                self.binding_data: np.expand_dims(binds, -1)
                }

        predictions = self.sess.run(self.predict_op, feed_dict=feed_dict)
        roc_auc = roc_auc_score(binds, predictions) 
        return binds, predictions, roc_auc 
    
    def train(self, data_directory):
        """
        Fully train the model given the input parameters. Runs for the input / default number of epochs. 
        Mostly just wraps the optimize_step function. 

        args:
            - data_directory: the input data directory. The data is split using the input metadata.  
        """
        train_data, test_data = utils.split_data(data_directory, self.validation_chunk, include_all=self.include_all_data, excluded_chunks=self.exclude_indices, n_folds=self.n_folds, seed=self.data_split_seed)
        data_iterator = utils.training_iterator(train_data, chunk_size=self.chunk_size)

        self.costs = list() 
        for epoch in tqdm(range(self.epochs), desc="Epochs"):
            epoch_costs = self.optimize_step(next(data_iterator))
            self.costs.append(epoch_costs)
            
            if (epoch + 1) % self.validate_step == 0:
                _, _, auc = self.predict(test_data)
                print("AUC: {0} (Epoch: {1})".format(epoch+1, auc))

        self.final_predictions = self.predict(test_data)

    def save_model(self, output_directory):
        """
        Save model outputs from a trained model. The output files are of the form: <id>.<spec>.<ext>
        The id is defined as the current time (as a string). The spec defines the output, either metadata
        or the domain type of output. ext is the file extension, either a numpy file or a json file.

        args:
            - output_directory: a string defining a directory to output the saved model to.  
        """
        import time        
        
        model_output_id = str(time.time()).replace("-", "_")
        
        # Metadata. Saved as dict (to JSON file) with four keys defining model type, 
        # model format, the results of the last prediction, and parameters
        results = {
                "binds": list(map(int, self.final_predictions[0])),
                "predictions": list(map(float, self.final_predictions[1])),
                "auc": float(self.final_predictions[2])
                }
        output = {
                "Model Type": "HSM / ID",
                "Model Format": list(self.model_spec), 
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
        model_file = os.path.join(output_directory, "{0}.{1}.{2}.npz".format(model_output_id, self.model_spec.domain_type, self.model_spec.peptide_type))
        dw, pw, iw, bias = self.sess.run([self.domain_weights, self.peptide_weights, self.interaction_weights, self.bias]) 
        np.savez(model_file, domain_weights=dw, peptide_weights=pw, interact_weights=iw, bias=bias)
