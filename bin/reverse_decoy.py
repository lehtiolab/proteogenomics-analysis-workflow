#!/usr/bin/env python3
from Bio import SeqIO
import sys


with open(sys.argv[1]) as fp, open('decoy_{}'.format(sys.argv[1]), 'w') as wfp:
  for target in SeqIO.parse(fp, 'fasta'):
    decoy = target[::-1] 
    decoy.description = decoy.description.replace('ENS', 'decoy_ENS')
    decoy.id = 'decoy_{}'.format(decoy.id)
    SeqIO.write(decoy, wfp, 'fasta')
