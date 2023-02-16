# PepInt predictions

PepInt models were used to make predictions SH2 and PDZ for comparison to HSM. See `driver.sh` for commands executed. To use the bash script, the following variables need to be set:

* `domain`: Specify PBD family to make predictions for. Must be one of `SH2` and `PDZ`.
* `data_path`: Specify path to data to make predictions for. This must be raw, unaligned data provided in csv format for a single PBD and must contain the columns (with header):
```
Domain UniProt ID,Domain Sequence,Peptidic Sequence,Bound
```
Raw HSM data (`/data/data_without_processed_duplicates/raw_data/`) can be downloaded from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529).
* `pepint_path`: Specify location of compiled PepInt models. Download PDZPepInt from [http://www.bioinf.uni-freiburg.de/Software/PDZPepInt/PDZPepInt.tar.gz](http://www.bioinf.uni-freiburg.de/Software/PDZPepInt/PDZPepInt.tar.gz) and compile per instructions. Download SH2PepInt from [http://www.bioinf.uni-freiburg.de/Software/SH2PepInt/SH2PepInt.tar.gz](http://www.bioinf.uni-freiburg.de/Software/SH2PepInt/SH2PepInt.tar.gz) and compile per instructions (after `make clean`). 

In order to make predictions with the HSM data, each domain was mapped to its PepInt model, and these mappings are provided in `map_domain_pepint_model`. Either UniProt IDs (if only one domain of a specificied PBD family is found in the protein) or raw domain sequences are used for the mapping. For futher information about the PDZPepInt models, see [Kundu, et al. Cluster based prediction of PDZ-peptide interactions (2014)](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-S1-S5). For further information about the SH2PepInt models, see [Kundu, et al. Semi-Supervised Prediction of SH2-Peptide Interactions from Imbalanced High-Throughput Data (2013)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0062732).