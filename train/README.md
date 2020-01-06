# Description

This directory may be used to train / re-train new models using the HSM framework.  

# Usage

Code used for training (or re-training) models is located in the `train/` directory in this repository. The package should primarily be accessed via the script `train.py`. 

```bash
python train.py [OPTIONS] 
```
Additional options for using `train.py` may be listed using the `-h/--help` flag. 

The basic steps for training a new model are:
0. Pre-process domain-peptide interaction data.

By default, the training code assumes that pre-processed data are located at `data/train/`, which can be downloaded (see [Data section](#data)). New data must be passed explicitly to the code (see the next section). 

A script for converting domain-peptide interaction data into the appropriate format for use with the model is available at `convert_binding_data.py`. The input format for this script is a csv file (no header) with the format: 
```
Domain-Type,Aligned-Domain-Sequence,Peptide-Type,Aligned-Peptidic-Sequence
```
An example of the input file type is included with the downloaded data (`domain_peptide_train.csv`). Domain and peptide protein identifiers are typically UniProtKB IDs; however, this is not required. Domain type refers to the class of binding domain (e.g. SH2 domains). 


To process the data, run the command:
```bash
python convert_binding_data.py [INPUT DATA FILE] [OUTPUT DATA DIRECTORY]
```
The processed data is output to the `OUTPUT DATA DIRECTORY` argument with the data split into directories by domain-type. Each directory contains data in two formats: `tf-records/` and `hdf5-records`. Additional options are detailed using the `-h/--help` flag. 

1. Train new models

Typically, models should be trained using the command:

```bash
python train.py [VALIDATION INDEX] (-d [DOMAIN ...] | -a) (--generic | --shared_basis | --hierarchical) 
```

The `VALIDATION INDEX` denotes data chunk that is excluded from the training process. The next argument, `-d [DOMAIN ...] | -a`, identifies the domains used in training the model. `-d` (or `--domains`) specifies a single or a subset of domains to train on. `a` (or `--all_domains`) specifies use all domains available. The final argument `(--generic | --shared_basis | --hierarchical)` specifies the model type: `--generic` specifies HSM/ID , `--hierarchical` denotes HSM/D, `--shared_basis` denotes HSM/D models trained for a single domain.

Additional command-line options facilitate model training / optimization (*e.g.* regularization parameters) and are detailed with the help command. 

2. Predict and assess performance

Data for the training process are typically output to `train/outputs/`. Processing the output directory can be accomplished using the `assess_performance.py` script:

```bash
python assess_performance.py [INPUT DIRECTORY]
```
where `INPUT DIRECTORY` denotes the path to the previously output directory. To control model training output, it can be helpful to re-direct outputs using the `--output_directory` command when running `train.py`.

3. Finalize model

To predict a combined model (*i.e.* using all training data) add the `--include_all` flag to the code 

```bash
python train.py [TEST INDEX] --include_all 
```

Models (for use with the `predict` code) may be output using the `output_models.py` script:
```bash
python output_models.py [RESULTS FILE] [OUTPUT DIRECTORY]
```
where `RESULTS FILE` denotes a file output by the `train.py` script and `OUTPUT DIRECTORY` the directory to place the processed models into (each domain type occupies one model file).  

