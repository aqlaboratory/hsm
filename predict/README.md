# Description

This directory may be used to make novel domain-peptide or protein-protein interaction predictions using pre-trained models.

# Usage

Predictions are run through one of two scripts, `predict_domains.py` and `predict_proteins.py`, for either domain-peptide or protein-protein predictions. All domain-peptide interactions that comprise a protein-protein prediction are output from the `predict_proteins.py` script. By default, data is output to an `outputs/` directory; please create this directory before running. 

Before using the pre-trained models (released in `models/`, please refer to the **Pre-trained Models** section, below.

## Domain-Peptide Interaction Predictions

Code used for predicting domain-peptide interactions is located in the predict/ directory in this repository. The functionality should primarily be accessed via the `predict_domains.py` script.

```python
python predict_domains.py [INPUT DOMAINS METADATA] [INPUT PEPTIDES METADATA] [OPTIONS] 
```
Additional options for using either script may be listed using the `-h/--help` flag. 

The basic steps for predicting a new interaction is:
### 0. Pre-process data and models.

By default, the code assumes that models are located at `predict/models/`. Output model files should be the same as formatted by `output_models.py` in the `train/` directory.

Pre-processed input domain and peptide metadata is available from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529), specifically in `data/ppi_data/metadata`.

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
python predict_proteins.py [-ppi [INPUT PPI PAIRS]] [OPTIONS] 
```
Additional options for using either script may be listed using the `-h/--help` flag. 

## 0. Pre-process data and models.

By default, the `predict_proteins.py` script also assumes models are located at `predict/models/`. The same model files may be used in both domain-peptide and protein-protein interaction prediction. To use new models, the same steps to specify the new models must be passed to `predict_proteins.py` as described above.

In addition, the models *require* metadata files that describe the domain and peptide composition of proteins. By defualt, the code assumes that pre-processed metadata for domains and peptidic sites identified in the human proteome is available at `../data/ppi_data/metadata/`. Metadata (`data/ppi_data/metadata`) can be downloaded from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529). New metadata must be passed explicitly to the code (see the next section).

## 1. Run predictions

Predictions can be computed using the described script:

```
python predict_proteins.py [--ppi_pairs [INPUT PPI PAIRS]] [OPTIONS] 
```
The `INPUT PPI PAIRS` option (passed using `--ppi_pairs`) passed to the code denotes a csv file containing the proteins to predict. These pairs should be formatted as a csv file where each line contains a pair of protein IDs (`<ID 1>,<ID 2>`). These IDs should reference IDs in the metadata files. If no pairs are passed, all valid pairs are returned. Different metadata files may be passed in using the `--domain_metadata` and `--peptide_metadata` options.  

## Pre-trained Models
Released with this codebase are pre-trained models for using HSM/D to make novel predictions. Pre-trained model files can be verified with the following MD5 hashes:

```
MD5 (amino_acid_ordering.txt) = c6458d6bd0cd9c5660632557b12cb28a
MD5 (hsm_pretrained/Kinase_TK.npz) = 746c6f4901ec85965faca7fe47ede548
MD5 (hsm_pretrained/PDZ.npz) = 0e34d0b3894c28f077122d82a2ffc910
MD5 (hsm_pretrained/PTB.npz) = d1c0e700e5f4f395eed38683433802a3
MD5 (hsm_pretrained/PTP.npz) = a8a1cf42dac2cfd16ec8757d20743ea7
MD5 (hsm_pretrained/SH2.npz) = 650303cb9088ec31a8305e73e9420b11
MD5 (hsm_pretrained/SH3.npz) = 802cb0a9fccb1fcfce9480993407a98a
MD5 (hsm_pretrained/WH1.npz) = 60e72ea081ed8c955ed5ec88fdf4b7b7
MD5 (hsm_pretrained/WW.npz) = 6aaf0473dba3aa32fe1336332503c09b
MD5 (hsm_pretrained/model_formats.csv) = 4e6f1842a6307ed5891ad9a9e9c67aa0
```

The code assumes that the input sequences are correctly aligned. Models were trained using the domain alignments released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `data/ppi_data/metadata/domain_metadata.csv`. MSAs constructed for each PBD family to obtain these alignments are also released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `results/data/alignment/`.

Peptides are aligned as follows:

| **Peptidic-type** | **Alignment** |
| ---------- | ---------------- |
| phosphosite | A 15 residue window aligned on the central phosphotyrosine. For example, `AAAAAAAyAAAAAAA`.|
| C-terminal | A 8 residue window aligned (to the right) on the C-terminus. For example, `AAAAAAAA`. Note, alignment offsets should be to the right (*e.g.* `AAAAAA` would become `--AAAAAA`).|
| polyproline | Any residue sequence. polyproline peptides are "scanning", or the likelihood is computed over the whole sequence.|

`A` denotes any amino acid. 

Note, using a different domain alignment and / or peptidic alignment (for 'fixed' peptides) results in predictions that are meaningless within the context of the model.
