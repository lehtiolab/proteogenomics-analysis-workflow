#!/usr/bin/env python3

import sys
import getopt
import os
from collections import OrderedDict

class Peptide(object):
    def __init__(self,chr=None,strand=None,start=0,end=0,psm=0,content=None):
        self.start=start
        self.end=end
        self.strand=strand
        self.chr=chr
        self.content=content


distance  = 10000 # default distance to group peptides is 10kb

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python group_novpepToLoci.py --input novpep.table.txt --output novpep.table.Loci.txt --distance 10kb")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=','output=','distance='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--output': output_file=arg
        elif opt == '--distance':
            try:
                distance=int(arg.replace("kb",""))*1000
            except ValueError:
                print ("incorrect input for distance");sys.exit()
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

input1=open(input_file,"r")
output=open(output_file,"w")

peplist=[]

header = input1.readline().strip().split("\t")
header += ["Loci","Loci_info"]
output.write("\t".join(header)+"\n")

idx=header.index("chr")

for line in input1: # 4rd to 7th column are chr, start, end, strand
    row=line.strip().split("\t")
    chr=row[idx].replace("chr","").replace("X","23").replace("Y","24")
    chr = chr.replace('NA', '25')
    peplist.append(Peptide(chr=chr,start=int(row[idx+1]),end=int(row[idx+2]),strand=row[idx+3],content=row))

## sort by chr and then start cor
peplist.sort(key=lambda x:list(map(int,(x.chr,x.start))))

print("total number of peptides",len(peplist))

loci_group=OrderedDict()
loci_num=1

for i in range(0,len(peplist)-1):
    if i==0:
        loci_group[loci_num]=[peplist[i]]
    
    if min(abs(peplist[i].start-peplist[i+1].end), abs(peplist[i].end-peplist[i+1].start))<distance:
        loci_group[loci_num].append(peplist[i+1])
    else:
        loci_num+=1
        loci_group[loci_num]=[peplist[i+1]]

print("total number of coding loci after groupping",len(loci_group))

for loci in loci_group:
    peptide_list=loci_group[loci]
    loci_info = "supported by %d peptide" % len(peptide_list)
    
    for pep in peptide_list:
        row = pep.content + [str(loci),loci_info]
        output.write("\t".join(row)+"\n")

input1.close()
output.close()
