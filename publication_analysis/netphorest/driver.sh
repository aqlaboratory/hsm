#!/bin/bash

domain= # name of PBD to make predictions for
data_path= # path to location of raw (unaligned) data

NETPHOR= # netphorest 2.1 executable

if [ ! -d $domain ]; then mkdir $domain; fi
cd $domain
   # write fasta for netphorest predictions
   # using raw data file since has UniProt IDs & netphorest independently aligns
   data=$data_path/$domain.csv
   python ../write_fasta.py $data

   # make predictions with netphorest
   if [ "$domain" == "Kinase_TK" ]; then
      classifier='KIN'
   else
      classifier=$domain
   fi
   cat peptides.fasta | $NETPHOR | grep $classifier > netphorest_predictions.tab
   
   # process predictions
   python ../process_netphorest_predictions.py netphorest_predictions.tab $data ../map_domain_netphorest_model/${domain}_netphorest_model.csv
cd ..
