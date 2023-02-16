"""
Process NetPhorest predictions

created by: Julia R Rogers
created on: 7-15-22

"""

import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

def create_parser():
   parser = argparse.ArgumentParser()
   parser.add_argument("netphorest_predictions", help="Output from NetPhorest.")
   parser.add_argument("data", help="Dataset in csv format for which predictions were made.")
   parser.add_argument("domain_model_map", help="Mapping between domain UniProt IDs in dataset and NetPhorest model.")

   return parser

def load_domain_model_map(csv):
   ds = np.loadtxt(csv, skiprows=1, delimiter=',', dtype=str)
   model_map = {d[0]: d[1] for d in ds if d[2] == '0'}
   seq_map = {d[-1]: d[1] for d in ds if d[2] == '1'}
   model_map.update(seq_map)
   return model_map, ds[:,0].flatten(), ds[:,1].flatten()

def load_netphorest_predictions(fname, model_list):
   colnames = ['Name', 'Position', 'Residue', 'Peptide',
               'Method', 'Tree', 'Classifier', 'Model', 'Posterior', 'Prior']
   cols = ['Name', 'Position', 'Model', 'Posterior']
   df = pd.read_csv(fname, sep='\t', names=colnames, header=None, usecols=cols)
   return df[df['Model'].isin(model_list)]

def match_predictions_dataset(dfpred, df, model_map):
   preds = []
   remove_rows = []
   for i, row in df.iterrows():
      if row['Domain UniProt ID'] in model_map.keys():
         model = model_map[row['Domain UniProt ID']]
      else: model = model_map[row['Domain Sequence']]
      name = row['Peptidic Sequence']
      pos = name.index('y')+1
      m = dfpred[(dfpred['Name'] == name) & (dfpred['Position'] == pos)
                 & (dfpred['Model'] == model)]
      if len(m) == 0:
         print("Netphorest failed to make a prediction for", row['Domain UniProt ID'], name, pos, model)
         remove_rows.append(i)
      else:
         preds.append(m.iloc[0]['Posterior'])
   df = df.drop(index=remove_rows)
   df['NetPhorest Posterior'] = preds
   return df

if __name__=='__main__':
   parser = create_parser()
   args = parser.parse_args()

   # mapping between individual domains and NetPhorest models in dataset
   model_map, domain_list, model_list = load_domain_model_map(args.domain_model_map)

   # load NetPhorest predictions
   dfpred = load_netphorest_predictions(args.netphorest_predictions, model_list)

   # load dataset
   df = pd.read_csv(args.data, header=0)
   df = df[df['Domain UniProt ID'].isin(domain_list)]

   # match NetPhorest predictions to dataset
   df = match_predictions_dataset(dfpred, df, model_map)
   basename = os.path.splitext(os.path.basename(args.data))[0]
   df.to_csv('%s_netphorest_predictions.csv' % basename, index=False)

   # construct ROC and calc AUC
   ys = df['Bound']
   preds = df['NetPhorest Posterior']
   auc = roc_auc_score(ys, preds)
   fprs, tprs, _ = roc_curve(ys, preds)
   with open('%s_roc.csv' % basename, 'w') as f:
      f.write('FPR,TPR,AUROC: %f\n' % auc)
   pd.DataFrame({'FPR': fprs, 'TPR': tprs}).to_csv('%s_roc.csv' % basename, index=False, header=False, mode='a')
   


