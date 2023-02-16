# PSSM models

HSM was compared to PSSMs constructed for each PBD family. To construct each PSSM model, the raw sequences of each domain in a PBD family modeled by HSM were first clustered based on intersequence distance defined by the PAM120 substitution matrix and computed after pairwise alignemnt. Then, the peptide sequences, aligned in the same manner as for HSM, were used to construct the PSSM for each cluster.

For a PBD family, PSSM model construction, evaluation using k-fold cross-validation (including tuning of clustering threshold) is performed with

```
python pssm.py [PREPROCESSED_DATA] [DOMAIN_SEQ_MAP] [PBD_FAMILY] [OPTIONS]
```
where

* `[PREPROCESSED_DATA]` is the data used for PSSM construction and evaluation. Data should be provided in a csv file with the columns (but no header):
```
Domain-Type,Aligned-Domain-Sequence,Peptide-Type,Aligned-Peptidic-Sequence,Bound
```
i.e. Data provided must be in the same format as the HSM pre-processed data. Pre-processed HSM data (`/data/data_without_processed_duplicates/preprocessed_raw_data.csv`) can be downloaded from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529).
* `[DOMAIN_SEQ_MAP]` is a mapping between raw, unaligned domain sequences and aligned, pre-processed sequences used for HSM provided in csv format. These are provided for the PBDs modeled in HSM in `map_domain_raw_preprocessed/`.
* `[PBD_FAMILY]` is the PBD family to model. One of `PTB`, `SH2`, `Kinase_TK`, `WW`, `WH1`, `PTP`, `SH3`, and `PDZ`.

Additional options (number of folds and number of processes to use) can be specified with the appropriate flags listed with the `-h/--help` flag.

A bash script for constructing PSSMs for all PBD families modeled by HSM is also provided (`driver.sh`).

## Requirements

* Python 3.9
* NumPy 1.23
* SciPy 1.8
* Pandas 1.4.2
* BioPython 1.79
* Scikit-Learn 1.0.2

Note that not all versions of BioPython include the PAM120 subsitution matrix by default. If needed, download from [ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/](ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/) and add to `<path_to_biopython>/Align/substitution_matrices/data/` where `<path_to_biopython>` is the path to the BioPython package being used.