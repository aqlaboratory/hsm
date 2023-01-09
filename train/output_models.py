"""
Reformats trained models for making new predictions

"""

import os
import numpy as np
import argparse

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("trained_model_directory", type=str, help="Directory with trained models to reformat.")
    parser.add_argument("model_specs", type=str, help="Model specifications contained with vectorized data used to train models.")
    parser.add_argument("-o", "--output_directory", type=str, default='hsm_pretrained', help='Directory to output reformatted models. Default=hsm_pretrained')
    return parser

def write_model_formats(args):
   # create model_specifications.csv
   domain_fixed = {'Kinase_TK': 1, 'WW': 0, 'PDZ': 1, 'WH1': 0, 'SH3': 0,
                   'PTB': 1, 'SH2': 1, 'PTP': 1}
   training_specs = np.loadtxt(args.model_specs, delimiter=',', usecols=range(4), dtype=str)
   with open(os.path.join(args.output_directory, 'model_formats.csv'), 'w') as f:
      for row in training_specs:
         dom = row[0]
         f.write('%s,%s,%s,%s,%i,%s.npz\n' % (*row, domain_fixed[dom], dom))

def reformat_model_params(args):
   # output model params with updated file names
   # and 'interact_weights' renamed to 'interaction_weights'
   npzs = [s for s in os.listdir(args.trained_model_directory) if s.endswith('npz')]
   for npz in npzs:
      d = dict(np.load(os.path.join(args.trained_model_directory, npz)))
      d['interaction_weights'] = d.pop('interact_weights')
      fname = os.path.join(args.output_directory, '%s.npz' % npz.split('.')[2])
      np.savez(fname, **d)


if __name__=='__main__':
   parser = create_parser()
   args = parser.parse_args()

   # create output directory if doesn't exist
   if not os.path.exists(args.output_directory):
      os.makedirs(args.output_directory)
   else:
      if len(os.listdir(args.output_directory)) > 0:
         overwrite = input("\n%s is not empty!!\nContinue anyways and overwrite contents? [y/n]\n" % args.output_directory)
         if overwrite != 'y' and overwrite != 'Y':
            sys.exit()

   write_model_formats(args)
   reformat_model_params(args)

