#!/usr/bin/env python3

import sys
import re
import getopt


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print ("Warning! wrong command, please read the mannual in Readme.txt.")
    print ("Example: python parse_spectrumAI_out.py --spectrumAI_out specAI_file --input input_psms.txt --output output_filename")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['spectrumAI_out=',
                                                         'input=',
                                                         'output='])
    for opt, arg in options:
        if opt == '--spectrumAI_out': input1_file=arg
        elif opt == '--input':input2_file=arg
        elif opt == '--output': output_file=arg
        else:
            print(("Warning! Command-line argument: %s not recognized. Exiting...") % opt); sys.exit()

input1=open(input1_file,"r")  ## SpectrumAI output
input2=open(input2_file,"r")  ## novelpep table, peptide sequence is in first column
output=open(output_file,"w")
header2=input2.readline().strip().split("\t")
header2 += ["SpectrumAI_result"]
output.write("\t".join(header2)+"\n")

header1=input1.readline().split("\t")
index1=header1.index("Peptide")
try:
    index2=header1.index("flanking_ions_support")
    index3=header1.index("status")
except ValueError:
    print ("the SpecturmAI output is empty")
    sys.exit()

specAI_result={}  # found with b,y ion support

for line in input1:
    row=line.strip().split("\t")
    pep=re.sub("[\W\d]","",row[index1].strip())
    try:
        if row[index3]=="checked":
            specAI_result[pep]=row[index2]
    except IndexError:
        print(("the line doesn't have the right number of columns"), line)


n1=0
n2=0
for line in input2: # peptide sequence is in first column
    row=line.strip().split("\t")
    pep=re.sub("[\W\d]","",row[0].strip())
    n1+=1
    if pep in specAI_result:
        row.append(specAI_result[pep])
        if specAI_result[pep]=="YES":
            n2+=1
    else:
        row.append("NA")

    output.write("\t".join(row)+"\n")

input1.close()
input2.close()
output.close()

print(("%d out of %d single substitution novel peptides passed SpectrumAI curation" % (n2,n1)))
