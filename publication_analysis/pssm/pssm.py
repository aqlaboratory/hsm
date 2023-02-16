"""
PSSM model

"""

import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.cluster.hierarchy import linkage, fcluster
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo, MultipleSeqAlignment, substitution_matrices
from Bio import pairwise2
from multiprocessing import Pool


def create_parser():
   parser = argparse.ArgumentParser()
   parser.add_argument("data", help="Dataset in csv format.")
   parser.add_argument("domain_seq_mapping", help="Mapping between aligned domain sequences and raw sequences in csv format.")
   parser.add_argument("domain_type", help="Domain-Type to model.")
   parser.add_argument("--nfolds", type=int, default=8, help="Number of folds for cross-validation. Default=8")
   parser.add_argument("--nprocs", type=int, default=1, help="Number of processes to use. Default=1")

   return parser

class PSSM():
   def __init__(self, domain_seqs, procs=1):
      self.seqs = np.array(domain_seqs)
      self.nseqs = len(domain_seqs)

      # construct distance matrix with distance defined by PAM120 subst matrix
      self.pam120_matrix = substitution_matrices.load('pam120')
      mat = np.zeros((self.nseqs, self.nseqs))
      p = Pool(processes=procs)
      # construct diagonal with scores first
      scores = p.starmap(self.align_score, zip(range(self.nseqs), range(self.nseqs)))
      np.fill_diagonal(mat, scores)
      # calc off-diagonal distances
      ijs = [(i, j) for i in range(self.nseqs-1) for j in range(i+1,self.nseqs)]
      scores = p.starmap(self.align_score, ijs)
      for score, ij in zip(scores, ijs):
         i, j = ij
         max_score = max(mat[i,i], mat[j,j])
         if max_score == 0: mat[i,j] = mat[j,i] = 1.0 
         else: mat[i,j] = mat[j,i] = 1.0 - (score * 1.0 / max_score)
      np.fill_diagonal(mat, 0)

      # range of sequence distances
      self.dist_range = (np.min(mat[mat > 0]), np.max(mat))
      # construct tree using hierarchical clustering
      self.Z = linkage(mat[np.triu_indices(mat.shape[0], k=1)], method='average')

   def align_score(self, i, j):
      return pairwise2.align.globaldx(self.seqs[i], self.seqs[j], self.pam120_matrix, score_only=True)

   def build_pssm(self, seqs):
      aln = MultipleSeqAlignment([SeqRecord(Seq(s)) for s in seqs])
      aln_info = AlignInfo.SummaryInfo(aln)
      pssm = aln_info.pos_specific_score_matrix().pssm
      for x, row in pssm:
         norm = sum([row[k] for k in row.keys()])
         for k in row.keys():
            try: row[k] = row[k]/norm
            except: row[k] = 0.0
      return pssm

   def train(self, data_df, threshold):
      # assign seqs to clusters
      c = fcluster(self.Z, threshold, criterion='distance')
      nclusters = max(c)
      cluster_domains = {i: self.seqs[np.where(c == i)] for i in range(1,nclusters+1)}
      self.map_to_cluster = {self.seqs[i]: c[i]-1 for i in range(self.nseqs)}
      
      # build PSSM from binding data for each cluster
      self.pssms = []
      for i in range(1, nclusters+1):
         df = data_df[data_df['Domain Sequence'].isin(cluster_domains[i])]
         if len(df) == 0:
            self.pssms.append(None)
            continue
         seqs = df['Aligned-Peptidic-Sequence'].tolist()
         self.pssms.append(self.build_pssm(seqs))

   def pssm_likelihood(self, seq, c):      
      if self.pssms[c] == None:
         return 0
      else:
         p = [self.pssms[c][i][1][aa] for i, aa in enumerate(seq) if aa in self.pssms[c][i][1].keys()]
         return np.prod(p)

   def predict(self, data_df):
      preds = []
      doms = data_df['Domain Sequence'].tolist()
      peps = data_df['Aligned-Peptidic-Sequence'].tolist()

      for dom, pep in zip(doms, peps):
         c = self.map_to_cluster[dom]
         preds.append(self.pssm_likelihood(pep, c))

      return preds

def load_domain_seq_map(csv):
   # Aligned-Domain-Sequence,Domain Sequence
   ds = np.loadtxt(csv, skiprows=1, delimiter=',', dtype=str, usecols=(1,2))
   return {d[0]: d[1] for d in ds}

def split_data(df, validation_chunk, n_folds=8, seed=0):
   np.random.seed(seed)
   randomized = np.random.permutation(len(df))
   chunks = np.array_split(randomized, n_folds)
   vndxs = chunks[validation_chunk]
   tndxs = [i for i in range(len(df)) if i not in vndxs]

   return df.iloc[tndxs].copy(), df.iloc[vndxs].copy()

def tune(train_df, val_df, model):
   bound_df = train_df[train_df['Bound'] == 1]
   ys = val_df['Bound']

   opt_auc = 0
   opt_pred = None
   opt_threshold = 1
   step = 0.01
   start = int(model.dist_range[0]/step)*step
   end = int(model.dist_range[1]/step+1)*step
   #print(model.dist_range, start, end)
   print("Trying thresholds in range [%10g, %10g)" % (start, end))
   for i in np.arange(start, end, step):
      model.train(bound_df, i)
      preds = model.predict(val_df)
      auc = roc_auc_score(ys, preds)
      print("Threshold:", i, "AUC:", auc)
      if auc > opt_auc:
         opt_auc = auc
         opt_pred = preds
         opt_threshold = i

   print("Optimal threshold: %8g with auc: %8.4f\n" % (opt_threshold, opt_auc))
   val_df['PSSM likelihood'] = opt_pred
   return val_df

if __name__=='__main__':
   parser = create_parser()
   args = parser.parse_args()

   colnames = ['Domain-Type','Aligned-Domain-Sequence',
               'Peptide-Type','Aligned-Peptidic-Sequence','Bound']
   df = pd.read_csv(args.data, header=None, names=colnames)
   df = df[df['Domain-Type'] == args.domain_type]
   print("PSSM model for %s" % args.domain_type)

   # initialize model with domain sequences
   dom_seq_map = load_domain_seq_map(args.domain_seq_mapping)
   model = PSSM(list(dom_seq_map.values()), args.nprocs)
   # map aligned domain sequences to (raw) domain sequences
   df['Domain Sequence'] = df['Aligned-Domain-Sequence'].map(dom_seq_map)

   # cross-validation with tuning of clustering threshold
   val_dfs = []
   for i in range(args.nfolds):
      print("\nTuning/training fold %i" % i)
      train_df, val_df = split_data(df, i, n_folds=args.nfolds)
      val_dfs.append(tune(train_df, val_df, model))

   val_df = pd.concat(val_dfs)
   val_df.to_csv('%s_pssm_predictions.csv' % args.domain_type, index=False)

   # construct ROC and calc AUC
   ys = val_df['Bound']
   preds = val_df['PSSM likelihood']
   auc = roc_auc_score(ys, preds)
   fprs, tprs, _ = roc_curve(ys, preds)
   with open('%s_roc.csv' % args.domain_type, 'w') as f:
      f.write('FPR,TPR,AUROC: %f\n' % auc)
   pd.DataFrame({'FPR': fprs, 'TPR': tprs}).to_csv('%s_roc.csv' % args.domain_type, index=False, header=False, mode='a')
   


