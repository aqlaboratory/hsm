# Description

This directory may be used to train / re-train new models using the HSM framework.  

**Note**, by default, results are output to an `output/` directory by default. This should be created.   

# Usage

Code used for training (or re-training) models is located in the `train/` directory in this repository. The package should primarily be accessed via the script `train.py`. 

```bash
python train.py [OPTIONS] 
```
Additional options for using `train.py` may be listed using the `-h/--help` flag. 

The basic steps for training a new model are:
## 0. Pre-process domain-peptide interaction data.

A script for converting domain-peptide interaction data into the appropriate format for use with the model is available at `convert_binding_data.py`. The input format for this script is a csv file (no header) with the format: 
```
Domain-Type,Aligned-Domain-Sequence,Peptide-Type,Aligned-Peptidic-Sequence
```
The input data used to train all HSM models is available on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) (`data/data_without_processed_duplicates/preprocessed_raw_data.csv`). Domain type refers to the class of binding domain (*e.g.* SH2). Peptide type refers to the class of peptide (*e.g.* phosphosite)


To process the data, run the command:
```bash
python convert_binding_data.py [INPUT DATA FILE] [OUTPUT DATA DIRECTORY]
```
The processed data is output to the `OUTPUT DATA DIRECTORY` with "*vectorized*" data split into directories of the form `<domain_type>_<peptide_type>`. Vectorized refers to converting domain and peptide sequences into one-hot representations.

By default, the training code expects processed data to be found at `../data/data_without_processed_duplicates/vectorized_data`.

Additional options are detailed using the `-h/--help` flag. 

## 1. Train new models

From the command line, models should be trained using the command:

```bash
python train.py [VALIDATION INDEX] (-d [DOMAIN ...] | -a) (--hsm_id | --hsm_d) 
```

The `VALIDATION INDEX` denotes data chunk that is excluded from the training process. The next argument, `-d [DOMAIN ...] | -a`, identifies the domains used in training the model. `-d` (or `--domains`) specifies a single or a subset of domains to train on. `a` (or `--all_domains`) specifies use all domains available. The final argument `(--hsm_id | --hsm_d)` specifies the model type: `--hsm_id` specifies HSM/ID and `--hsm_d` denotes HSM/D. Note, multiple domains may only be passed to the `--hsm_d` option. Both the `--hsm_id` and `--hsm_d` options may accept single-domain arguments. 

Additional command-line options facilitate model training / optimization (*e.g.* regularization parameters) and are detailed with the help command.

Models may also be imported from the `models/` package. For examples of how to configure these models, please reference the `train.py` script. 

### Data splits.
To simplify the model's implementation, this code implements data splitting by setting a single random seed (`np.random.seed()`). Consequently, to train over a single data split, pass the same random seed (using the `--chunk_seed` option) to the model.



## 2. Utilize outputs

Data and results for the training process are typically output to the specified output directory (defaults to `outputs/`). Raw models are output as NumPy array archives (.npz, output using `np.savez()`). The arrays have keys:
```
domain_weights, peptide_weights, interact_weights, basis_weights, bias
``` 
that describe the model. 

The output metadata completely describes the parameters of the model. It is output as a JSON file containing four keys:
```
Model Type, Model Formats, Results, Parameters
```
These keys have the model IDs (`Model Type, Model Formats`) as well as input parameters (`Parameters`) and the previous performance and predictions, output in the results file. 


## 3. Finalize models

Typically, model training will involve two steps - training with different folds withheld and training on the combined data set. 

To test transferability, it is typical to identify a test fold (*e.g.* 5) and then optimize hyperparameters relative to a second fold (*e.g.* 2). To do this, the user could start by excluding folds 2 and 5 and test on fold 2:
```bash
python train.py 2 --exclude_chunks 5 (--hsm_id | --hsm_d) --domains [DOMAINS...] [OPTIONS]
```
In this step, folds 2 and 5 are excluded from the data and performance is assessed on fold 2. Having optimized hyperparameters, the model might be re-trained using all folds by fold 5:
```bash
python train.py 5 (--hsm_id | --hsm_d) --domains [DOMAINS...] [OPTIONS]
```

By default, the validation fold is excluded from training. Sometimes, it is necessary to include this fold (e.g. for a final model). To do so, use the `--include_all_data` argument:
```bash
python train.py 5 --include_all_data (--hsm_id | --hsm_d) --domains [DOMAINS...] [OPTIONS]
```
The model will still be assessed against fold 5. Note, this means that the training process 'sees' the data in fold 5 before that assessment; therefore, this is not an appropriate way to measure performance. 


# Training of published HSM models

Specifics of hyperparameter tuning and training used for the published HSM models is described below. `train_foldtune_grid.py` was used to tune HSM/ID hyperparameters using grid search and train a model on each fold split. `train_foldtune_rand.py` was used to tune HSM/D hyperparameters using random search and train a model on each fold split. Detailed descriptions and results from hyperparameter tuning are provided in `Analyze_training.ipynb`.

## HSM/ID models

For each PBD family, an HSM/ID model was trained and evaluated using 8-fold cross-validation. For each fold, hyperparameters were selected based on a grid search. Each fold was trained using the command:

```bash
mkdir fold[i]
python train_foldtune_grid.py [i] --hsm_id -d [DOMAIN] -o fold[i]
```
where `i` is the index of the test fold and `DOMAIN` is the PBD family modeled. Other options available to specify grid search parameters are detailed using the `-h/--help` flag. Defaults were used for all published models.

## HSM/D models

An HSM/D model was trained and evaluated using 8-fold cross-validation. Two rounds of random hyperparameter search were performed. In the first round, models were trained using the command:

```bash
mkdir fold[i]
python train_foldtune_rand.py [i] --hsm_d -a --min_log_learning_rate -5 --max_log_learning_rate -3 --lower_lambda_diff -3 --upper_lambda_diff 1.5 --ntrials 50 --epochs 120 --validate_step 5 --lambda_params "SH2:1e-05 Kinase_TK:8.762500000000001e-06 PTP:0.0001 SH3:4.3750000000000005e-06 WW:8.7625e-05 PDZ:4.3750000000000005e-06 WH1:2.665e-05 PTB:2.8e-05" -o fold[i]
```
where `i` is the index of the test fold.

In the second round, models were trained using the command:

```bash
mkdir fold[i]
python train_foldtune_rand.py [i] --hsm_d -a --min_log_learning_rate -4.5 --max_log_learning_rate -3 --lower_lambda_diff -0.7 --upper_lambda_diff 0.7 --ntrials 30 --epochs 200 --validate_step 5 --lambda_params "SH2:5.80156e-07 Kinase_TK:3.02538e-06 PTP:1.13193e-05 SH3:6.9499e-07 WW:2.09734e-06 PDZ:3.02872e-07 WH1:3.0459e-05 PTB:0.000169522" -o fold[i]
```
where `i` is the index of the test fold.

A final HSM/D model was trainined on all data for use in HSM/P using the command:

```bash
mkdir final
python train.py 0 --include_all_data --hsm_d -a -r 0.00030746435 -e 94 -v 94 --lambda_params "SH2:8.3553e-07 Kinase_TK:4.28949e-06 PTP:2.15278e-05 SH3:6.88512e-07 WW:2.7692e-06 PDZ:4.10312e-07 WH1:2.97909e-05 PTB:0.000178861" -o final
```

Output from the final HSM/D model training was reformatted for use in HSM/P using the command:

```bash
python output_models.py final <DATADIR>/vectorized_data/models_specification.csv
```
where `DATADIR` is the directory with the "*vectorized*" data used for training all models. `models_specification.csv` is generated when the raw data is "*vectorized*" for training as described above in **0. Pre-process domain-peptide interaction data.**