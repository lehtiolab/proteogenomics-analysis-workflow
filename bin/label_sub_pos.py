#!/usr/bin/env python3

import sys
import re
import os
import getopt

################  Comand-line arguments ################
peptide_column="Peptide" #default name to look for peptide column in input file
protein_id_column = "Protein" #default name to look for protein accessions in input file
splitchar = "_" #default splitting by '_'

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python lab_sub_pos.py --input_psm PSM_filename --output output_filename")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input_psm=','output=','splitchar='])
    for opt, arg in options:
        if opt == '--input_psm': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--splitchar': splitchar=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

input1=open(input_file,"r") # vardb search psm table

header1=input1.readline().strip().split("\t")

try:
    pep_index=header1.index(peptide_column)
    pro_index=header1.index(protein_id_column)
except ValueError:
    print("Peptide, Protein columns are not in input table"); sys.exit()

output=open(output_file,"w") # 1 mismatch novel peptide PSM

newheader=header1+["sub_pos"]
output.write("\t".join(newheader)+"\n")

for line in input1:
    row=line.strip().split("\t")
    pro=row[pro_index]

    sub_pos="NA"
    acc=pro.split(";")[0]
    acc=re.sub(r'\([^)]*\)', '', acc) # remove text within parentheses
    
    if acc[:6]=="CanPro":
        sub_pos=acc.split("_")[-1]
    elif acc[:6]=="COSMIC":
        sub_pos=acc.split(":")[-1]
    else:
        splitheader = acc.split(splitchar)
        if len(splitheader) == 1:
            continue  # No split, no position
        try:
            int(splitheader[-1])  # something else found instead of position
        except:
            continue
        else:
            sub_pos = splitheader[-1]

    row.append(sub_pos)
    output.write("\t".join(row)+"\n")

input1.close()
output.close()
