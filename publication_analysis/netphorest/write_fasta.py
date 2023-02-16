"""
write fasta files of peptides for use with NetPhorest

created on: 7-15-22
created by: Julia Rogers
"""

import sys, os
import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1], header=0)

pep_seqs = df['Peptidic Sequence'].drop_duplicates().to_list()

with open('peptides.fasta', 'w') as f:
   for s in pep_seqs:
      f.write('>%s\n%s\n' % (s,s.upper()))

