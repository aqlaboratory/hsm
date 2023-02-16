import os, csv
import numpy as np
import json
from collections import *

def load_results(directory, domain, nfolds):
   """
   Load predictions and target values from CV metadata outputs
   Returns numpy array with labels in first column and predictions in second
   """
   data = [[], []]

   for fold in range(nfolds):
      folddir = os.path.join(directory, 'fold%i' % fold)
      fname = [s for s in os.listdir(folddir) if s.endswith('metadata.json')][0]
      d = json.load(open(os.path.join(folddir, fname)))
      # HSM/D
      if isinstance(d['results'], list):
         fname = [s for s in os.listdir(folddir) if s.endswith('mean.metadata.json')][0]
         print(domain, fname)
         d = json.load(open(os.path.join(folddir, fname)))
         for result in d['results']:
            if result[0][0] == domain:
               data[0].extend(result[1]['binds'])
               data[1].extend(result[1]['predictions'])
               print(fold, len(result[1]['predictions']))
               break
      # HSM/ID
      else:
         print(domain, fname)
         data[0].extend(d['results']['binds'])
         data[1].extend(d['results']['predictions'])
         print(fold, len(d['results']['predictions']))

   return np.array(data)

def load_hsmd(fname):
   with open(fname, 'r') as f:
      lines = f.readlines()
   # each row of data: P1_uid, P2_uid, list of domain-peptide predictions
   data = [json.loads(l) for l in lines]

   # list of tuples of domain model, predicted probability for each PPI 
   dpi_probs = [[(dpi[0][0], dpi[2]) for dpi in d[2]] for d in data]

   prots2d = np.array([[d[0], d[1]] for d in data], dtype=str)
   prots1d = np.array([','.join(sorted(prot)) for prot in prots2d])

   return prots1d, prots2d, dpi_probs, data

def load_model_threshold_fpr(directory):
   "Loads data from HSM/D ROCs for each PBD family needed to calculate p-vals."

   fs = [f for f in os.listdir(directory) if 'fpr_threshold.csv' in f]

   mdict = {}
   for f in fs:
      k = f.split('_')[0]
      if k == 'Kinase': k = 'Kinase_TK'
      mdict[k] = np.loadtxt(os.path.join(directory, f), skiprows=1, delimiter=',')
   return mdict

def load_prots(fname):
   "Load protein pairs from csv file."
   prots2d = np.loadtxt(fname, delimiter=',', usecols=(0,1), dtype=str)
   prots2d = strip_isoform_uniprotid(prots2d)
   prots1d = np.array([','.join(sorted(prot)) for prot in prots2d])

   return prots1d, prots2d

def strip_isoform_uniprotid(a):
   "Strips isoform information from UniProt IDs provided in numpy array"
   return np.vectorize(lambda x: x.split('-')[0])(a)

def load_protein_metadata(domain_metadata_fp='../data/ppi_data/metadata/domain_metadata.csv', peptide_metadata_fp='../data/ppi_data/metadata/peptide_metadata.csv'):
    "Loads HSM protein metadata"
    mdata = defaultdict(list)
    for uid, seq, dtype in csv.reader(open(domain_metadata_fp, 'r')):
        mdata[uid].append((dtype, seq))
    for uid, seq, ptype, _ in csv.reader(open(peptide_metadata_fp, 'r')):
        mdata[uid].append((ptype, seq))
        
    return mdata

def get_protein_composition(uid, metadata=None):
    if metadata is None: metadata = load_protein_metadata()
    return metadata[uid]

def get_ppi_composition(uid1, uid2, interactions, p, threshold=0.05, threshold_operator=np.less_equal, metadata=None):
   "Get PPI composition that only includes PBDs and SLiMs involved in at least one interaction given criteria."
   rc1 = get_protein_composition(uid1, metadata=metadata)
   rc2 = get_protein_composition(uid2, metadata=metadata)
    
   components = set()
   for ip, i in zip(p, interactions):
      if threshold_operator(ip, threshold):
         for j in range(2):
            components.add(tuple(i[j]))
   c1 = [c for c in rc1 if c in components]
   c2 = [c for c in rc2 if c in components]
    
   # return compositions of each prot, set of tuples (pbd/slim_type, seq)
   return c1, c2

def output_data_json(fname, data):
   "Outputs each entry in data list as json string."
   with open(fname, 'w+') as f:
      for d in data:
         f.write(json.dumps(d) + "\n")

