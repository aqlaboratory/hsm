# Description

This directory may be used to train / re-train new models using the HSM framework.  

**Note**, by default, data are assumed to be at `../data/` downloaded as detailed in the `README.md` file of the parent directory. In addition, results are output to an `output/` directory by default. This should be created.   

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
An example of the input file type is included with the downloaded data (`data/train/hsm_preprocessed/preprocessed_raw_data.csv`). Domain type refers to the class of binding domain (*e.g.* SH2 domains). Peptide type refers to the class of peptide (*e.g.* phosphosite)


To process the data, run the command:
```bash
python convert_binding_data.py [INPUT DATA FILE] [OUTPUT DATA DIRECTORY]
```
The processed data is output to the `OUTPUT DATA DIRECTORY` with "*vectorized*" data split into directories of the form `<domain_type>_<peptide_type>`. Vectorized refers to converting domain and peptide sequences into one-hot representations. 

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

Data and results for the training process are typically output to `train/outputs/`. Raw models are output as NumPy array archives (.npz, output using `np.savez()`). The arrays have keys:
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
