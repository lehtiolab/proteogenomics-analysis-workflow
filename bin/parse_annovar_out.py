#!/usr/bin/env python3

import sys
import getopt


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python prepare_annovar_input.py --input novpep.tab.txt --annovar_out novpep.variant_function --output novpep_anovar.txt")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','output=', 'annovar_out='])
    for opt, arg in options:
        if opt == '--annovar_out': annovar_file=arg
        elif opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()


input1=open(annovar_file,"r") # annovar output
input2=open(input_file,"r") # novel pep table
output=open(output_file,"w")

category={}
while True:
    line1=input1.readline()
    line2=input1.readline()
    if not line2:break
    row1=line1.strip().split("\t")
    row2=line2.strip().split("\t")
    pep=row1[-1].replace("Comments:Seq=","")
    if row1[0]==row2[0]: # if the Nterm and Cterm of the peptide match to same functional class
        category[pep]=[row1[0],row1[1]]
    else:
        fun=row1[0]+"-"+row2[0]
        category[pep]=[fun,row1[1]+";"+row2[1]]

print("%d peptides returned with anovar annotation" % len(category))

input2.readline()
header=["Peptide","anovar_category","associated_gene"]
output.write("\t".join(header)+"\n")
for line in input2:
    row=line.strip().split("\t")
    pep=row[0]
    
    anovar_cat = "NA"
    associated_gene="NA"
    if pep in category:
        anovar_cat=category[pep][0]
        associated_gene=category[pep][1]

    output.write("\t".join([pep,anovar_cat,associated_gene])+"\n")

input1.close()
input2.close()
output.close()






