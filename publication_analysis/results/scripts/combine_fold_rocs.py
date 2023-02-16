"""
Combine ROC curves from cross-fold validation

created by: Julia R Rogers
created on: 8-6-21

"""

import os, sys, argparse
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from utils import load_results

def create_parser():
   """
   Create argument parser
   """
   parser = argparse.ArgumentParser()
   parser.add_argument("training_directory", help="Directory with output from CV.")
   parser.add_argument("domain", help="Domain to construct ROC for.") 
   parser.add_argument("--nfolds", type=int, default=8, help="Number of folds for cross-validation. Default: 8")

   return parser


if __name__=='__main__':
   parser = create_parser()
   args = parser.parse_args()

   print(args.training_directory)
   data  = load_results(args.training_directory, args.domain, args.nfolds)

   auc = roc_auc_score(data[0,:], data[1,:])
   fprs, tprs, thresholds = roc_curve(data[0,:], data[1,:])

   with open('%s.csv' % args.domain, 'w') as f:
      f.write('FPR,TPR,AUROC: %f\n' % auc)
   pd.DataFrame({'FPR': fprs, 'TPR': tprs}).to_csv('%s.csv' % args.domain, index=False, header=False, mode='a')

   pd.DataFrame({'FPR': fprs[1:], 'Thresholds': thresholds[1:]}).to_csv('%s_fpr_threshold.csv' % args.domain, index=False)
   
