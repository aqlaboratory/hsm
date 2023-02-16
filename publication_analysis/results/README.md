# Description

This directory contains code and notebooks used to analyze HSM model results and create all figures in the publication. Data generated from (or needed for) analysis is available on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529), specifically in `results/data/`; upon download, the subdirectory `data/` should be placed in this directory to reproduce the analysis.

## HSM/ID and HSM/D ROCs and model calibration (Figure 2a, Supplementary Figure 3, Supplementary Table 3)

### Receiver operating characteristic (ROC) curves

Receiver operating characteristic (ROC) curves are calculated for each PBD family with `combine_fold_rocs.py` provided in `scripts/` using the command:

```
python scripts/combine_fold_rocs.py [TRAINING DIRECTORY] [PBD] [OPTIONS]
```
where `[TRAINING DIRECTORY]` is the directory with the output from cross-fold validation and `[PBD]` is the PBD family to construct the ROC curve for. Additional options can be listed with the flag `-h/--help`. HSM training data is released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `training/` and serves as an example of how the directory should be structured: `[TRAINING DIRECTORY]` should contain subdirectories named `fold%i`, which each contain output from training a single fold `i` trained as described in **Training of published HSM models** in the directory `train/` of this repo.

All calculated ROCs are released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `results/data/perf/hsmid/` and `results/data/perf/hsmd/`.

### Model calibration

Model calibration curves are calculated for each PBD family with `experimental_proportions.py` provided in `scripts/` using the command:

```
python scripts/experimental_proportions.py [TRAINING DIRECTORY] [PBD] [OPTIONS]
```
where arguments are as specified above for calculating ROC curves.

All calculated calibration curves are released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `results/data/perf/hsmid_exp_proportions/` and `results/data/perf/hsmd_exp_proportions/`.

### Figure 2a and Supplementary Figure 3

After calculating (or downloading the ROCs and calibration curves), Figure 2a and Supplementary Figure 3 are created with `Fig_2a_S3.ipynb`. To fully recreate these figures, ROCs for all external models and PSSMs must also be recomputed (see `publication_analysis/netphorest`, `publication_analysis/pepint`, and `publication_analysis/pssm` directories of this repo) or download from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `netphorest/`, `pepint/`, and `pssm/` and place in this directory's parent directory.

### Supplementary Table 3

Significance testing was performed using the DeLong test for all PBD models. p-values are calculated with `Table_S3_pvals_delong.ipynb` using a fast Python implementation provided on [https://github.com/yandexdataschool/roc_comparison](https://github.com/yandexdataschool/roc_comparison). To reproduce these calculations, download predictions for each model from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `training/`, `netphorest/`, `pepint/`, and `pssm/` in addition to raw training data in `data/data_without_processed_duplicates/raw_data/` and place in this direcotry's parent directory.

## HSM/P evaluations and analysis of proteome-scale predictions (Figures 2b & 3, Supplementary Figure 4, and Supplementary Table 5)

Evaluations, comparisons to high-throughput experiments, and analysis of newly predicted PPIs were made with `Fig_2b_Table_S5_inputs_Fig_3_S4_analysis_hsmp_proteome_predictions.ipynb`.

### Figure 2b and Supplementary Table 5

To reproduce this analysis, recreate Figure 2b and obtain numbers reported in Supplementary Table 5, download all HSM/P predictions and experimental datasets from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `predictions/` and `data/ppi_data/ppi_expt_comparisons/` and place in this directory's parent directory. To compute HSM/D and HSM/P p-values, the thresholds used to construct the ROCs must be computed (thresholds are saved when ROCs are computed as described above) or downloaded from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `results/data/perf/hsmid/` and `results/data/perf/hsmd/`.

Results from analysis of HSM/P predictions are released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `results/data/perf/hsmp/`. This includes the list of all HSM predicted PPIs at a FDR of 0.01 and notes which ones have been experimentally confirmed.

### Figure 3 and Supplementary Figure 4

After using `Fig_2b_Table_S5_inputs_Fig_3_S4_analysis_hsmp_proteome_predictions.ipynb` to identify interaction mechanisms of PPIs predicted by HSM/P at a FDR of 0.01 and confirmed by experiments conducted after or at the same time as HSM training data was reported, Figure 3 and Supplementary Figure 4 are created with Mathematica notebook `Fig_3_S4.nb` (which requires data formatted as done in `Fig_2b_Table_S5_inputs_Fig_3_S4_analysis_hsmp_proteome_predictions.ipynb`). To create these required inputs, download the preprocessed domain and peptide metadata used with HSM/P from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `data/ppi_data/metadata/` and place in this directory's parent directory.

## Insights into PBD–peptide interaction mechanisms from HSM/D (Figures 4 & 5 and Supplementary Figures 5 & 6)

Biophysical analysis of HSM/D models and creation of Figures 4 & 5 and Supplementary Figures 5 & 6 was done with `Fig_4_5_S5_S6.ipynb`. To reproduce this analysis and recreate the Figures, download model weights, domain alignments, HCK SH3 structure, and PyMOL session files from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529) in `results/data/weights/`, `results/data/alignment/`, and `results/data/structures/`.

## HSM/P derived PPI networks (Figure 6 and Supplementary Figure 7)

PPI networks are constructed with `Fig_6_S7.ipynb`. To recreate Figure 6 and Supplementary Figure 7, download all HSM/P predictions for the human proteome and annotated with p-values using `Fig_2b_Table_S5_inputs_Fig_3_S4_analysis_hsmp_proteome_predictions.ipynb` (`results/data/perf/hsmp/ppi_predictions_pvals.csv`) from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529). Additionally,  download the preprocessed domain and peptide metadata used with HSM/P from figshare in `data/ppi_data/metadata/` and place in this directory's parent directory.

# Requirements

Versions used for the publication are noted.

* Python 3.9
* NumPy 1.23
* StatsModels 0.13.2
* SciPy 1.8
* Matplotlib 3.5.1
* Seaborn 0.11.2
* Graph-Tool 2.45
* Scikit-Learn 1.0.2
* BioPython 1.79
* Pandas 1.4.2
* Mathematica 13.2
* PyMOL 2.5.3