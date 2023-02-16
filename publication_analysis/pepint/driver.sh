#!/bin/bash

domain= # name of PBD to make predictions for
data_path= # path to location of raw (unaligned) data

pepint_path= # path to directory with compiled PepInt models
PEPINT=$pepint_path/${domain}PepInt

if [ ! -d $domain ]; then mkdir $domain; fi
cd $domain
   cwd=`pwd`
   # using raw data file since has UniProt IDs & pepint independently aligns
   data=$data_path/$domain.csv
   
   # write fasta for pepint predictions
   python ../write_fasta.py $data

   # run PepInt
   cd $PEPINT
      cmd=$PEPINT/${domain}PepInt.sh
      stdout=`$cmd $cwd/peptides.fasta`
      name=`echo $stdout | cut -c 12-`
      mv $name $cwd/
   cd $cwd

   # process predictions
   python ../process_pepint_predictions.py $name $data ../map_domain_pepint_model/${domain}_pepint_model.csv
cd ..
