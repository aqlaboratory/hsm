domains=( PTB SH2 Kinase_TK WW WH1 PTP SH3 PDZ )

mapping_dir=map_domain_raw_preprocessed
data= # specify path to data, i.e. preprocessed_raw_data.csv

mkdir -p output
cd output

for domain in ${domains[@]}; do
   mfile=$mapping_dir/${domain}_matched_domseqs.csv
   python ../pssm.py ../$data ../$mfile $domain --nprocs 6
done

cd ..


