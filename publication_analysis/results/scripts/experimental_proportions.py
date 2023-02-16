"""
Parse model output for Fig S3 comparison to epxeriment

"""

import os, sys, argparse
import numpy as np
import pandas as pd
import csv

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

def segment_data(data, ofile, step_size=0.05):
    """
    Computes empirical proportion over data 

    Inputs:
        data: numpy array of binary values and predictions for each data pt
        step_size: step size to compute experimental proportion over. 
    """

    outputs = list()

    # Iterates over the range [0,1] splitting it into appropriate step sizes.
    # Note, the mean of the binary values is equivalent to the proportion of positives.
    a = np.linspace(0,1,(int(1/step_size) + 1))
    for i,j in zip(a, a[1:]):
        mask = np.logical_and(i <= data[1,:], data[1,:] < j)
        outputs.append((i, np.mean(data[0,:][mask])))

    with open(ofile, 'w+') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(["Interval Start", "Empirical Proportion"])
        writer.writerows(outputs)

if __name__=='__main__':
   parser = create_parser()
   args = parser.parse_args()

   data = load_results(args.training_directory, args.domain, args.nfolds)

   segment_data(data, '%s.csv' % args.domain) 

   
