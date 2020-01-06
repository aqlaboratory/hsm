# Description

This directory may be used to make novel domain-peptide or protein-protein interaction predictions using pre-trained models.

# Usage

Predictions are run through one of two scripts, `predict_domains.py` and `predict_proteins.py`, for either domain-peptide or protein-protein predictions. All domain-peptide interactions that comprise a protein-protein prediction are output from the `predict_proteins.py` script. By default, data is output to an `outputs/` directory; please create this directory before running. 

## Domain-Peptide Interaction Predictions

Code used for predicting domain-peptide interactions is located in the predict/ directory in this repository. The functionality should primarily be accessed via the `predict_domains.py` script.

```python
python predict_domains.py [INPUT DOMAINS METADATA] [INPUT PEPTIDES METADATA] [OPTIONS] 
```
Additional options for using either script may be listed using the `-h/--help` flag. 

The basic steps for predicting a new interaction is:
### 0. Pre-process data and models.

By default, the code assumes that models are located at `predict/models/` and pre-processed data, which can be downloaded from [figshare (doi:10.6084/m9.figshare.11520552)](https://figshare.com/articles/Pre-processed_data_-_Git_Repo_-_HSM/11520552), should be available at `data/predict`. New data must be passed explicitly to the code (see the next section). Output model files should be the same as formatted by `output_models.py` in the `train/` directory. 

Input domains files should have the format:
```
Domain-Protein-Identifier,Aligned-Domain-Sequence,Domain-Type
```
and input peptide files:

```
Peptide-Protein-Identifier,[Aligned-?]Peptide-Sequence,Peptide-Type,Is-Fixed
```

The peptide should be aligned for fixed models (e.g. phosphosites) and not-aligned for non-fixed models. 


### 1. Run predictions

Predictions can be computed using the described script:

```python
python predict_domains.py [INPUT DOMAINS METADATA] [INPUT PEPTIDES METADATA] [OPTIONS] 
```

Domain-peptide interactions are computed for all valid pairs (*e.g.* pairs that have an associated model). The two major options, `-m/--models` and `--model-format`, are useful when using newly trained models. `-m/--models` describes a directory listing new model files (each model file should specify the associated domain type). `--model-format` describes the model metadata as a csv file:
```
Domain-Type,Peptide-Type,Domain-Alignment-Length,Peptide-Alignment-Length,Peptide-Alignment-Is-Fixed,Model-Filename
```

The domain and peptide alignment lengths refer to the domain / peptide alignments used to train the model. The is-fixed option describes if the interaction should be computed using a scanning scope (*e.g.* for a poly-proline site in the paper) or a fixed scope (*e.g.* for a phosphosite or C-terminus site in the paper) using 0/1 for scanning / fixed. The model-filename lists the associated file name (note, not file path).

## Protein-Protein Interaction Predictions

Code used for predicting protein-protein interactions is located in the predict/ directory in this repository. The functionality should primarily be accessed via the `predict_proteins.py` script.

```python
python predict_proteins.py [-p [INPUT PPI PAIRS]] [OPTIONS] 
```
Additional options for using either script may be listed using the `-h/--help` flag. 

## 0. Pre-process data and models.

By default, the `predict_proteins.py` script also assumes models are located at `predict/models/` and pre-processed data, which can be downloaded via [figshare (doi:10.6084/m9.figshare.11520552)](https://figshare.com/articles/Pre-processed_data_-_Git_Repo_-_HSM/11520552), are available at `data/metadata`. New data must be passed explicitly to the code (see the next section). The same models files may be used in both domain-peptide and protein-protein interaction prediction. To use new models, the same steps to specify the new models must be passed to `predict_proteins.py`. In addition, the models requiire metadata files (by default, stored in `data/metadata`) that describe either the domain or peptide composition of proteins. Metadata are formatted as Python dictionaries (stored as pickle'd files) with the format: 

## 1. Run predictions

Predictions can be computed using the described script:

```python
python predict_proteins.py [--ppi_pairs [INPUT PPI PAIRS]] [OPTIONS] 
```
The `INPUT PPI PAIRS` option (passed using `--ppi_pairs`) passed to the code denotes a csv file containing the proteins to predict. These pairs should be formatted as a csv file where each line contains a pair of protein IDs (`<ID 1>,<ID 2>`). These IDs should reference IDs in the metadata files. If no pairs are passed, all valid pairs are returned. Different metadata files may be passed in using the `--domain_metadata` and `--peptide_metadata` options.  
