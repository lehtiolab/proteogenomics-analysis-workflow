#!/usr/bin/env python3
from Bio import SeqIO
import sys


with open(sys.argv[1]) as fp, open('td_concat.fa', 'w') as wfp:
  for target in SeqIO.parse(fp, 'fasta'):
    SeqIO.write(target, wfp, 'fasta')
    decoy = target[::-1] 
    decoy.description = decoy.description.replace('ENS', 'decoy_ENS')
    decoy.id = 'decoy_{}'.format(decoy.id)
    SeqIO.write(decoy, wfp, 'fasta')
