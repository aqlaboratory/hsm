"""
Process PepInt predictions

created by: Julia R Rogers
created on: 10-21-22

"""

import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

def create_parser():
   parser = argparse.ArgumentParser()
   parser.add_argument("pepint_predictions", help="Output from PepInt.")
   parser.add_argument("data", help="Dataset in csv format for which predictions were made.")
   parser.add_argument("domain_model_map", help="Mapping between domain UniProt IDs in dataset and PepInt model.")

   return parser

def load_domain_model_map(csv):
   ds = np.loadtxt(csv, skiprows=1, delimiter=',', dtype=str)
   model_map = {d[0]: d[1] for d in ds if d[2] == '0'}
   seq_map = {d[-1]: d[1] for d in ds if d[2] == '1'}
   model_map.update(seq_map)
   return model_map, ds[:,0].flatten(), ds[:,1].flatten()

def load_pepint_predictions(fname, model_list):
   # SH2 model makes predictions for all Y possible and note which in Position
   # peptide aligned with 2 upstream, 4 downstream
   remove_rows = []
   try:
      cols = ['Seq-ID:', 'Position:', 'Sequence:', *model_list]
      df = pd.read_csv(fname, sep='\t', usecols=cols)
      # remove any peptide sequences that were not properly aligned
      for i, row in df.iterrows():
         if len(row['Sequence:']) < 7 or row['Sequence:'][2] != 'Y': remove_rows.append(i)
   # SH3 and PDZ don't have position column
   except:
      cols = ['Seq-ID:', 'Sequence:', *model_list]
      df = pd.read_csv(fname, sep='\t', usecols=cols)

      # SH3 model appends considered range to Seq-ID
      if len(df.iloc[0]['Seq-ID:'].split('-')) == 3:
         df['Seq-ID:'] = df['Seq-ID:'].apply(lambda d: d.split('-')[0])
      else:
         # remove any peptide sequences that were not properly aligned
         for i, row in df.iterrows():
            if len(row['Sequence:']) < 5: remove_rows.append(i)
   if len(remove_rows) > 0:
      print("Removed %i data points since PepInt could not properly align" % len(remove_rows))
   return df.drop(remove_rows)

def match_predictions_dataset(dfpred, df, model_map):
   preds = []
   remove_rows = []
   seqs_wo_models = []
   for i, row in df.iterrows():
      # determine which PepInt model to use
      if row['Domain UniProt ID'] in model_map.keys():
         model = model_map[row['Domain UniProt ID']]
      else:
         try: model = model_map[row['Domain Sequence']]
         except:
            seqs_wo_models.append((row['Domain UniProt ID'], row['Domain Sequence']))
            remove_rows.append(i)
            continue
      # peptide sequence
      name = row['Peptidic Sequence']
      # if SH2, check which prediction to use based on pY location
      if 'Position:' in dfpred.columns:
         pos = name.index('y')+1
         posname = 'Y%i' % pos
         m = dfpred[(dfpred['Seq-ID:'] == name) & (dfpred['Position:'] == posname)]
      else:
         m = dfpred[dfpred['Seq-ID:'] == name]
      if len(m) == 0 or model not in m.columns:
         print("PepInt failed to make a predictions for", row['Domain UniProt ID'], name, model)
         remove_rows.append(i)
      else:
         preds.append(m[model].item())
   df = df.drop(index=remove_rows)
   df['PepInt Score'] = preds

   print("No model for %i dom seqs" % len(set(seqs_wo_models)))
   for entry in set(seqs_wo_models):
      print("No model for: %s, %s" % (entry[0], entry[1]))

   return df

if __name__=='__main__':
   parser = create_parser()
   args = parser.parse_args()

   # mapping between individual domains and PepInt models in dataset
   model_map, domain_list, model_list = load_domain_model_map(args.domain_model_map)

   # load PepInt predictions
   dfpred = load_pepint_predictions(args.pepint_predictions, model_list)

   # load dataset
   df = pd.read_csv(args.data, header=0)
   df = df[df['Domain UniProt ID'].isin(domain_list)]

   # match PepInt predictions to dataset
   df = match_predictions_dataset(dfpred, df, model_map)
   basename = os.path.splitext(os.path.basename(args.data))[0]
   df.to_csv('%s_pepint_predictions.csv' % basename, index=False)

   # construct ROC and calc AUC
   ys = df['Bound']
   preds = df['PepInt Score']
   auc = roc_auc_score(ys, preds)
   fprs, tprs, _ = roc_curve(ys, preds)
   with open('%s_roc.csv' % basename, 'w') as f:
      f.write('FPR,TPR,AUROC: %f\n' % auc)
   pd.DataFrame({'FPR': fprs, 'TPR': tprs}).to_csv('%s_roc.csv' % basename, index=False, header=False, mode='a')
   

   # AUC per domain
   wav = 0
   total = 0
   print("\n\nAUC's per domain")
   udoms = df['Domain Sequence'].unique()
   for ud in udoms:
      ddf = df.loc[df['Domain Sequence'] == ud]
      uid = ddf['Domain UniProt ID'].unique()[0]
      ys = ddf['Bound'].to_numpy()
      preds = ddf['PepInt Score']
      print("%s, %s" % (uid, ud))
      posfrac = sum(ys == 1)/len(ddf)
      try:
         auc = roc_auc_score(ys, preds)
         print("AUC: %8.4f, %i pts (%8.4f), pos frac %8.4f\n" % (auc, len(ddf), len(ddf)/len(df), posfrac))
         wav += len(ddf)*auc
         total += len(ddf)
      except:
         y = ys[0]
         auc = sum(preds == y)/len(ddf)
         print("Frac of predictions matching label %i: %8.4f, %i pts\n" % (auc, y, len(ddf)))
         
   print("Weighted av of domain AUC: %8.4f" % (wav/total))
   print("computed from %i data points" % total)

