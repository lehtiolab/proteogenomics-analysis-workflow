#!/usr/bin/env python3


'''
    the script is written by Mikael Hussius @ SciLifeLab, https://github.com/hussius/gff-phastcons-human.
    slightly modified to work after the bigwig file is downloaded locally.
'''
import pyBigWig as pw
import numpy as np
import sys


if len(sys.argv)<4:
    sys.exit("USAGE: python " + sys.argv[0] + "<GFF file with regions> <BigWig file> <output file>")


# Need to download hg19.100way.phastCons.bw first, use the command:
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw


infile = sys.argv[1]
bw_file = sys.argv[2]
outfile = sys.argv[3]

bw = pw.open(bw_file)
oF = open(outfile, "w")

header = ["Bare peptide","phastcon_max_score","phastcon_mean_score"]
oF.write("\t".join(header)+"\n")

pep_dic = {}

if not infile.endswith("gff") and not infile.endswith("gff3"):
    sys.exit("The region file needs to be a GFF(3) file!")

for line in open(infile):
    if not line.startswith("chr"):
        continue
    fields = line.strip().split()
    (chr, start, end, pept) = (fields[0], fields[3], fields[4], fields[8])
    if not pept.startswith("Parent="): continue
    pept = pept.replace("Parent=","") # Remove "Parent="
    
    try:
        values = bw.values(chr, int(start), int(end))
        max_val = np.max(values)
        mean_val = np.mean(values)
    except:
        print("Encountered error for line: " + line.strip())
        max_val = -1
        mean_val = -1
    
    if pept not in pep_dic:
        pep_dic[pept]=[max_val, mean_val]
    else:
        pep_dic[pept][0]=max(max_val,pep_dic[pept][0])
        pep_dic[pept][1]=np.mean([mean_val,pep_dic[pept][1]])

bw.close()

for pept in pep_dic:
    oF.write("%s\t%f\t%f\n" % (pept,pep_dic[pept][0],pep_dic[pept][1]))

oF.close()
