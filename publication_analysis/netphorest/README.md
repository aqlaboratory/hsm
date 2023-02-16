# NetPhorest predictions

NetPhorest 2.1 was used to make predictions PTB, SH2, PTP, and TKs for comparison to HSM. See `driver.sh` for commands executed. To use the bash script, the following variables need to be set:

* `domain`: Specify PBD family to make predictions for. Must be one of `PTB`, `SH2`, `Kinase_TK`, and `PTB`.
* `data_path`: Specify path to data to make predictions for. This must be raw, unaligned data provided in csv format for a single PBD and must contain the columns (with header):
```
Domain UniProt ID,Domain Sequence,Peptidic Sequence,Bound
```
Raw HSM data (`/data/data_without_processed_duplicates/raw_data/`) can be downloaded from [figshare (doi:10.6084/m9.figshare.22105529)](https://doi.org/10.6084/m9.figshare.22105529).
* `NETPHOR`: Specificy netphorest 2.1 executable. Download `NetPhorest_human_2.1.zip` from [http://netphorest.science/download.shtml](http://netphorest.science/download.shtml) and compile per instructions.

In order to make predictions with the HSM data, each domain was mapped to its NetPhorest model, and these mappings are provided in `map_domain_netphorest_model`. Either UniProt IDs (if only one domain of a specificied PBD family is found in the protein) or raw domain sequences are used for the mapping. For futher information about the NetPhorest models, see [Miller, et al. Linear Motif Atlas for Phosphorylation-Dependent Signaling (2008)](https://www.science.org/doi/10.1126/scisignal.1159433) and [Horn, et al. KinomeXplorer: an integrated platform for kinome biology studies (2014)](https://www.nature.com/articles/nmeth.2968).