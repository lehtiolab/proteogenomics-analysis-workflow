#!/usr/bin/env python3

import sys
import os
import getopt
from Bio import SeqIO
import re

if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python parse_blastp_output.py --input input_filename --blastp_result peptide.blastp.out.txt --fasta Homo_sapiens.GRCh37.75.pep.all.fa --output output_filename")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                         'blastp_result=',
                                                         'fasta=',
                                                         'output='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--fasta': fasta_file=arg
        elif opt == '--blastp_result':blastp_file=arg
        elif opt == '--output': output_file=arg
        else:
            print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit()

input1= SeqIO.parse(fasta_file,"fasta") # blastp reference database
seqdb={}
for record in input1:
    seq=str(record.seq)
    if record.id not in seqdb:
        seqdb[record.id]=seq

print(len(seqdb))

input2= open(blastp_file,"r") # blastp output
input3= open(input_file,"r") # novel peptide tab txt
output= open(output_file,"w")

blastout={}
hits_dic={}
for line in input2:
    row=line.strip().split("\t")
    qid=row[0]
    sid=row[1]
    sseq=seqdb[sid]
    ident=row[2]
    peplen=int(row[3])
    alignlen=int(row[6])-int(row[5])+1
    sstart=int(row[7])
    send=int(row[8])
    mismatch=row[9]
    gap=row[12]
    alignseq=row[-3]
    evalue=float(row[-2])
    category="NA"
    single_sub_pos="NA"
    if sstart>3:
        Nterm_seq=sseq[sstart-4:sstart+2] #check up 3 amino acid before N-term of this peptide
    else:
        Nterm_seq=sseq[:sstart]

    if len(sseq)-send<3:
        Cterm_seq=sseq[send-1:]
    else:
        Cterm_seq=sseq[send-3:send+3]

    if alignlen==peplen:
        if float(ident)==100:
            category="match to known protein"
        
        elif int(gap)==0 and int(mismatch)==1:
            category="map to known protein with 1 aa mismatch"
            for i in range(peplen):
                if qid[i]!=alignseq[i]:
                    single_sub_pos=str(i+1)

        elif int(gap)==1 and int(mismatch)==0:
            category="map to known protein with 1 aa insertion"
        else:
            category="novelpep (map to known protein with more than 2 mismatched aa)"
    elif peplen-alignlen==1 and float(ident)==100:
        category="map to known protein with 1 aa deletion"

    else:
        category="novelpep (map to known protein with more than 2 mismatched aa)"
    
    if qid not in hits_dic:
        hits_dic[qid]=evalue
        blastout[qid]=[category,sid,ident,peplen,single_sub_pos,Nterm_seq,alignseq,Cterm_seq,alignlen,mismatch,gap]
    else:
        if evalue<hits_dic[qid]:
            hits_dic[qid]=evalue
            blastout[qid]=[category,sid,ident,peplen,single_sub_pos,Nterm_seq,alignseq,Cterm_seq,alignlen,mismatch,gap]

header=input3.readline().strip().split("\t")

header=header+["blastp_category","blastp_match","identity","peplen","sub_pos","Nterm-seq(3aa)","aligned_seq","Cterm-seq(3aa)","alignlen","mismatch","gap"]
output.write("\t".join(header)+"\n")

for line in input3:
    row=line.strip().split("\t")
    peptide=row[0]
    if peptide in blastout:
        results=blastout[row[0]]
        newrow=row+results
        output.write("\t".join(map(str,newrow))+"\n")
    else:
        newrow=row+["novelpep (no match to known protein found)","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"]
        output.write("\t".join(map(str,newrow))+"\n")


input1.close()
input2.close()
input3.close()
output.close()
