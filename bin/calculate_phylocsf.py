#!/usr/bin/env python3

'''
    the script is modified from Mikael Hussius @ SciLifeLab, https://github.com/hussius/gff-phylocsf-human
    
    download the following bigwig files first
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+0.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+1.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF+2.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-0.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-1.bw
    # wget https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg19/latest/PhyloCSF-2.bw
'''

import sys
import os
import pyBigWig as pw
import numpy as np

def predict_coding(vec):
    coding = "OTHER"
    for v in vec:
        if not v: continue
        if v > 0: coding = "CODING"
    return(coding)

if len(sys.argv)<4:
    sys.exit("USAGE: python " + sys.argv[0] + "<GFF file> <BigWig file path> <output file>")

infile = sys.argv[1]
bw_file_path = sys.argv[2]
outfile = sys.argv[3]

regs = []
chrom={}
starts={}
ends={}
peptide={}

for line in open(infile):
    if not line.startswith("chr"):
        continue
    fields = line.strip().split()
    (chr, start, end, pept) = (fields[0], fields[3], fields[4], fields[8])
    if not pept.startswith("Parent="): continue
    name = chr+":"+start+"-"+end
    chrom[name]=chr
    starts[name]=int(start)
    ends[name]=int(end)
    peptide[name]=pept.split("=")[1]
    regs.append(name)

scores = {}

rpathbase = os.path.join(bw_file_path,"PhyloCSF")

for rf in ["+0","+1","+2","+3","-0","-1","-2","-3"]:
    rpath = rpathbase + rf + ".bw"
    if os.path.isfile(rpath):
        sys.stderr.write("Searching PhyloCSF reading frame " + rf + "\n")
        bw = pw.open(rpath)
        frame_score = {}
        count = 0
        for r in regs:
            count += 1
            if(count % 50 ==0): sys.stderr.write('\tProcessed ' + str(count) + " peptides out of " + str(len(regs)) + "\n")
            sys.stderr.flush()
            try:
                score = bw.stats(chrom[r], starts[r], ends[r])[0]
            except RuntimeError:
                pass
            frame_score[r] = score
            scores[rf] = frame_score
        bw.close()
    else:
        sys.stderr.write("%s doesn't exist \n" % rpath)



output = open(outfile,"w")
output.write("\t".join(["Bare peptide","PhyloCSF+0.score","PhyloCSF+1.score","PhyloCSF+2.score","PhyloCSF-0.score","PhyloCSF-1.score","PhyloCSF-2.score","PhyloCSF_prediction"])+"\n")

pep_scores={}

for r in regs:
    scoreList = [scores["+0"][r], scores["+1"][r], scores["+2"][r], scores["-0"][r], scores["-1"][r], scores["-2"][r]]
    seq = peptide[r]
    if seq not in pep_scores:
        pep_scores[seq]=scoreList
    else: # this is to consider splice junction peptides which have two regions separated in gff file, we take mean phylocsf score of two regions
        for i in range(0,len(scoreList)):
            value = scoreList[i]
            if value is None and pep_scores[seq][i] is None:
                continue
            elif None in [value, pep_scores[seq][i]]:
                pep_scores[seq][i] = value if value else pep_scores[seq][i]
            else:
                pep_scores[seq][i] = (pep_scores[seq][i] + value)/2

for seq in pep_scores:
    scoreList = pep_scores[seq]
    row = [seq]+['NA' if x is None else str(x) for x in scoreList] + [predict_coding(scoreList)]
    output.write('\t'.join(row) + '\n')
