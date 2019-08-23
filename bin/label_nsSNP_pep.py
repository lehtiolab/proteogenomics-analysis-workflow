#!/usr/bin/env python3

import sys
import os
import getopt
from Bio import SeqIO


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print ("Warning! wrong command, please read the mannual in Readme.txt.")
    print ("Example: python lab_sub_pos.py --input pep.tab.txt --nsSNPdb nsSNP_pep.fa --output output_filename")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','nsSNPdb=','output='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--nsSNPdb': snp_db=arg
        elif opt == '--output': output_file=arg
        else:
            print ("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

handle= SeqIO.parse(snp_db,"fasta")
db_dic={}

for record in handle:
    if "rs" in record.description: # check if it has rs id in entry's description
        db_dic[record.description] = str(record.seq)

input1=open(input_file,"r") # novpep tab table

header= input1.readline().strip().split("\t")
header += ["from.SNPdb"]

output=open(output_file,"w") # output table
output.write('\t'.join(header) + '\n')

print ("searching SNPdb to see if any novel peptides derived from nsSNPs")
for line in input1:
    row=line.strip().split("\t")
    pep = row[0]
    
    query_output = "No match in SNP-DB"
    for pro in db_dic:
        if pep in db_dic[pro]:
            query_output=pro
            break

    row.append(query_output)
    output.write("\t".join(row)+"\n")

input1.close()
output.close()
handle.close()

