# hsm - Biophysical prediction of protein-peptide interactions and signaling networks using machine learning. 

<img align="left" src="misc/symbol_name.png" style="width: 25%; height: 25%"/> 

This repository implements the hierarchical statistical mechanical (HSM) model described in the paper [Biophysical prediction of protein-peptide interactions and signaling networks using machine learning.](https://doi.org/10.1038/s41592-019-0687-1) 

An **associated website** is available at [proteinpeptide.io](https://proteinpeptide.io). The website is built to facilitate interactions with results from the model including: (1) specific domain-peptide and protein-protein predictions, (2) the resulting networks, and (3) structures colored using the inferred energy functions from the model. Code for the website is available via the parallel repo: [aqlaboratory/hsm-web](https://github.com/aqlaboratory/hsm-web). Note that the results on the website were obtained using an [old model](#model-updates).

This file documents how this package might be [used](#usage), the [location of associated data](#data), and [other metadata](#reference). 

## Usage

The model was implemented in Python (>= 3.5) primarily using TensorFlow (>= 1.14) ([Software Requirements](#requirements)). To work with this repository, either download pre-processed data (see below) or include new data. The folder contains three major directories: `train/`, `predict/`, and `publication_analysis/`. Each directory is accompanied by a `README.md` file detailing usage. 

To train / re-train new models, use the `train.py` script in `train/`. To make predictions using a model, use one of two scripts, `predict_domains.py` and `predict_proteins.py`, for predicting either domain-peptide interactions or protein-protein interactions. Scripts are designed with a CLI and should be used from the command line: 

```bash
python [SCRIPT] [OPTIONS]
```

Options for any script may be listed using the `-h/--help` flag.

To reproduce analysis and figures presented in the paper [Biophysical prediction of protein-peptide interactions and signaling networks using machine learning.](https://doi.org/10.1038/s41592-019-0687-1), use the scripts in `publication_analysis/`.

Pre-trained models are released with this repo. An alternative use case would be to train / re-train a new model in the `train/` code and make new predictions using the `predict/` code. 

### Model updates

We identified an issue in the original datasets used to train the model published in [Biophysical prediction of protein-peptide interactions and signaling networks using machine learning.](https://doi.org/10.1038/s41592-019-0687-1). As of February 15, 2023, we have corrected the datasets (released on [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529)), and replaced the original models released with this repo with corrected ones. Please verify that you use the corrected models for all predictions (see documentation in `predict/`).

## Data

All associated data may be downloaded from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529).


## Requirements
- Python (>= 3.5)
- TensorFlow (1.14)
- numpy (1.18)
- scipy (1.4)
- scikit-learn (0.20)
- tqdm (4.41) (Progressbar. Not strictly necessary for functionality; needed to ensure package runs.)


## Reference
Please reference the associated publication:

Cunningham, J.M., Koytiger, G., Sorger, P.K., & AlQuraishi, M. "Biophysical prediction of protein-peptide interactions and signaling networks using machine learning." *Nature Methods* (2020). [doi:10.1038/s41592-019-0687-1](https://doi.org/10.1038/s41592-019-0687-1). ([citation.bib](misc/citation.bib))

See also, a **website** at [proteinpeptide.io](https://proteinpeptide.io) for exploring the associated analyses (code: [aqlaboratory/hsm-web](https://github.com/aqlaboratory/hsm-web)). 

## Funding

This work was supported by the following sources:

| **Funder** | **Grant number** |
| ---------- | ---------------- |
| NIH | U54-CA225088 |
| NIH | P50-GM107618 |
| DARPA / DOD | W911NF-14-1-0397 |

## License
This repository is released under an [MIT License](LICENSE)
