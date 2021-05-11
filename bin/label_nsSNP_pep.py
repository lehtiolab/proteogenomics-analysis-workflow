#!/usr/bin/env python3

import sys
import os
import getopt
import sqlite3
from Bio import SeqIO



if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print ("Warning! wrong command, please read the mannual in Readme.txt.")
    print ("Example: python lab_sub_pos.py --input pep.tab.txt --nsSNPdb nsSNP_pep.fa --output output_filename")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','nsSNPdb=','output=', 'dbfile=', 'minlen='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--nsSNPdb': snp_db=arg
        elif opt == '--dbfile': dbfn=arg
        elif opt == '--output': output_file=arg
        elif opt == '--minlen': minpeplen=arg
        else:
            print ("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

handle = SeqIO.parse(snp_db,"fasta")
whole_proteins = {}

for record in handle:
    if "rs" in record.description: # check if it has rs id in entry's description
        whole_proteins[record.description] = str(record.seq)
whole_proteins = {prot: seq.replace('L', 'I') for prot, seq in whole_proteins.items()}
 

input1=open(input_file,"r") # novpep tab table

header= input1.readline().strip().split("\t")
header += ["from.SNPdb"]

output=open(output_file,"w") # output table
output.write('\t'.join(header) + '\n')

print ("searching SNPdb to see if any novel peptides derived from nsSNPs")
conn = sqlite3.Connection(dbfn)

for line in input1:
    row=line.strip().split("\t")
    pep = row[0].replace('L', 'I')
    query_output = "No match in SNP-DB"
    cur = conn.execute('SELECT protid, pos FROM protein_peptides WHERE seq="{}"'.format(pep[:int(minpeplen)]))
    snpprots = [(protid, pos) for protid, pos in cur]
    for prot_id, pos in snpprots:
        protseq = whole_proteins[prot_id]
        if pep in protseq:
            query_output = prot_id
            break
    row.append(query_output)
    output.write("\t".join(row)+"\n")


input1.close()
output.close()
handle.close()

