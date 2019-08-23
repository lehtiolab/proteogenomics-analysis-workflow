#!/usr/bin/env python3

import sys
import re

input1 = open(sys.argv[1],"r") # blastp parse out novpep table
input2 = open(sys.argv[2],"r") # novel psm table

output = open(sys.argv[3],"w")

novpep_1mismatch = {}
cols = input1.readline().split("\t")
idx1 = cols.index("blastp_category")
idx2 = cols.index("sub_pos")

for line in input1:
    row = line.strip().split("\t")
    cat = row[idx1]
    pep = row[0]
    sub_pos = row[idx2]
    if cat=="map to known protein with 1 aa mismatch":
        novpep_1mismatch[pep] = sub_pos

print ("%d peptides map to known protein with 1 aa mismatch" % len(novpep_1mismatch))

header = input2.readline().strip().split("\t")
header += ["sub_pos"]

output.write("\t".join(header)+"\n")

idx = header.index('Peptide')
for line in input2:
    row = line.strip().split("\t")
    seq = re.sub('[^a-zA-Z]','',row[idx])
    if seq in novpep_1mismatch:
        row.append(novpep_1mismatch[seq])
        output.write("\t".join(row)+"\n")

input1.close()
input2.close()
output.close()
