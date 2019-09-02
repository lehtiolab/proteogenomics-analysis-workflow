#!/usr/bin/env python3

'''
    this script will map peptides to genome and output an peptide gff3 file with coordinates.
    which genome it maps to depends on the GTF input file provided.
    written by Yafeng Zhu @ Karolinska Institutet.  Email: yafeng.zhu@ki.se
'''

import sys
import os
import getopt
import numpy as np
import re
from Bio import SeqIO

class EXON(object):
    def __init__(self,number=0,gene=None,variant=None,chr=None,strand=None,start=0,end=0,length=0,trans_start=0,trans_end=0):
        self.gene=gene
        self.variant=variant
        self.number=number
        self.start=start  #chromosome start coordinate
        self.end=end  #chromosome end coordinate
        self.strand=strand
        self.chr=chr
        self.trans_start=trans_start
        self.trans_end=trans_end
    def length(self):
        self.length=self.trans_end-self.trans_start+1


def cal_trans_pos(exon_list): # calculate transcript position of exon start & end, exon_list is a list of exon objects
    strand=exon_list[0].strand
    if strand=="+":
        new_exon_list=sorted(exon_list,key=lambda x:x.start)
    else:
        new_exon_list=sorted(exon_list,key=lambda x:x.start,reverse=True)
    
    sumExonlength=0
    for exon in new_exon_list:
        exon_length=exon.end-exon.start+1
        exon.trans_start=1+sumExonlength
        exon.trans_end=exon.trans_start+exon_length-1
        sumExonlength+=exon_length
        
    return new_exon_list

def get_pep_cor(exon_object_list,n1,n2): # return peptide's chromosome start and end cor given peptide's trans_start (n1) and trans_end (n2)
    pep_chr=""
    pep_strand=""
    pep_chr_start=0
    pep_chr_end=0
    pep_start_exon=0
    pep_end_exon=0
    for i in range(len(exon_object_list)):
        exon=exon_object_list[i]
        if n1<=exon.trans_end and n1>=exon.trans_start:
            pep_chr=exon.chr
            pep_strand=exon.strand
            pep_start_exon=i+1
            if pep_strand=='+':
                pep_chr_start=exon.start+(n1-exon.trans_start)
            else:
                pep_chr_end=exon.end-(n1-exon.trans_start)

        if n2<=exon.trans_end and n2>=exon.trans_start:
            pep_chr=exon.chr
            pep_strand=exon.strand
            pep_end_exon=i+1
            if pep_strand=='+':
                pep_chr_end=exon.start+(n2-exon.trans_start)
            else: # chr_cor of n2 is pep_chr_start
                pep_chr_start=exon.end-(n2-exon.trans_start)

    return pep_chr,pep_strand,pep_chr_start,pep_chr_end,pep_start_exon,pep_end_exon

def parse_gtf(infile):
    dic={}
    with open(infile,"r") as infile_object:
        for line in infile_object:
            if line[0]!="#": # skip lines commmented out
                row=line.strip().split("\t")
                if row[2]=="exon":
                    attri_list=row[8].split(";")
                    transID=""
                    exon=EXON(start=int(row[3]),end=int(row[4]),chr=row[0],strand=row[6])
                    for attri in attri_list:
                        if "transcript_id " in attri:
                            transID=attri.strip().replace("transcript_id ","").replace('\"',"")
                        elif "transcript_id=" in attri:
                            transID=attri.strip().replace("transcript_id=","").replace('\"',"")
                
                    if transID not in dic:
                        dic[transID]=[exon]
                    else:
                        dic[transID].append(exon)
    return dic

def get_id(s):
    acclist=s.split(";")
    acc = acclist[0]
    prot_ID=acc.split("(")[0]
    trans_ID = ""
    if "PGOHUM_" in prot_ID:
        trans_ID=prot_ID.split("_")[1]
    elif "PGOHUM00000" in prot_ID:
        trans_ID=prot_ID.split("_")[0]
    elif "lnc-" in prot_ID:
        trans_ID=prot_ID.split("_")[0]
    elif "lncRNA" in prot_ID:
        trans_ID=prot_ID.split("_")[1]

    return prot_ID,trans_ID

################  Comand-line arguments ################


if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
    print("Warning! wrong command, please read the mannual in Readme.txt.")
    print("Example: python map_novelpeptide2genome.py --input input_filename --gtf VarDB.gtf --fastadb VarDB.fasta --tab_out file0.txt --fasta_out file1.fasta --gff3_out file2.gff3 --bed_out file3.bed")
else:
    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
                                                         'gtf=',
                                                         'fastadb=',
                                                         'tab_out=',
                                                         'fasta_out=',
                                                         'gff3_out=',
                                                         'bed_out='])
    for opt, arg in options:
        if opt == '--input': input_file=arg
        elif opt == '--gtf':gtf_file=arg
        elif opt == '--fastadb': db_file=arg
        elif opt == '--tab_out': tab_file=arg
        elif opt == '--fasta_out': fasta_file=arg
        elif opt == '--gff3_out': gff_file=arg
        elif opt == '--bed_out': bed_file=arg
        else:
            print(("Warning! Command-line argument: %s not recognized. Exiting..." % opt)); sys.exit()

print("reading GTF input file")
feature_dic=parse_gtf(gtf_file)
print(("number of unique transcripts in GTF file",len(feature_dic)))

seqdb = SeqIO.parse(db_file,'fasta')
seq_dic={}
for record in seqdb:
    if record.id[:6] in ["COSMIC","CanPro"]:
        continue
    else:
        seq_dic[record.id]=str(record.seq)

print(("number of lncRNA & pseudogene sequences in fasta file",len(seq_dic)))

input=open(input_file,'r') # peptide table with two columns, peptide sequence in first column, protein ID in second column

tab_output=open(tab_file,'w')
fasta_output=open(fasta_file,'w')
gff_output=open(gff_file,'w')
bed_output=open(bed_file,'w')

header=input.readline().split("\t")

pep_col = header.index("Peptide")
pro_col = header.index("Protein")

newheader=["Bare peptide","Peptide","Protein","chr","start","end","strand"]
tab_output.write("\t".join(newheader) + '\n')

non_mapped_pep=0
novpep_dic={}

for line in input:
    row=line.strip().split("\t")
    modpep = row[pep_col].strip()
    peptide=re.sub("[\W\d]","", modpep).upper()
    proteins = row[pro_col]
    
    if peptide not in novpep_dic:
        novpep_dic[peptide] = 1
        fasta_output.write(">%s\n%s\n" % (peptide,peptide))
    else:
        continue;
    
    if "chr" in proteins:
        # Only report first protein match coordinates, remove multiple-hits by BLAT later
        protein = proteins.split(';')[0]
        cor_list=protein.split("_")
        print(protein, cor_list)
        pep_chr=cor_list[0]
        pep_chr_start=cor_list[1]
        pep_chr_end=cor_list[2]
        pep_strand=cor_list[3]
        
        gff_format_line1=[pep_chr,"MS","mRNA",pep_chr_start,pep_chr_end,".",pep_strand,".","ID="+peptide]
        gff_format_line2=[pep_chr,"MS","CDS",pep_chr_start,pep_chr_end,".",pep_strand,"0","Parent="+peptide]
        gff_output.write("\t".join(map(str,gff_format_line1))+"\n")
        gff_output.write("\t".join(map(str,gff_format_line2))+"\n")
        
        newrow=[peptide,modpep,proteins,pep_chr,pep_chr_start,pep_chr_end,pep_strand]
        tab_output.write("\t".join(newrow)+"\n")
    
        continue;

    protein_id,transcript_id=get_id(proteins)
    frame=protein_id[-1]
    
    try:
        exons=feature_dic[transcript_id]
    except KeyError:
        non_mapped_pep+=1
        print(("KeyError",transcript_id,"doesn't exit in GTF input file"))
        continue;
    
    aa_seq=seq_dic[protein_id]
    pep_index=aa_seq.index(peptide)
    
    pep_trans_start=3*pep_index+1
    pep_trans_end=pep_trans_start+3*len(peptide)-1
    
    if frame=="2":
        pep_trans_start+=1
        pep_trans_end+=1
    elif frame=="3":
        pep_trans_start+=2
        pep_trans_end+=2

    exons=cal_trans_pos(exons)
    
    #print pep_trans_start,pep_trans_end
    pep_chr,pep_strand,pep_chr_start,pep_chr_end,pep_start_exon,pep_end_exon=get_pep_cor(exons,pep_trans_start,pep_trans_end)

    newrow=[peptide,modpep,proteins]+list(map(str,[pep_chr,pep_chr_start,pep_chr_end,pep_strand]))
    tab_output.write("\t".join(newrow)+"\n")

    bed_output.write("%s\t%s\t%s\tA\t-\tComments:Seq=%s\n" % (pep_chr,pep_chr_start,pep_chr_start,peptide))
    bed_output.write("%s\t%s\t%s\tA\t-\tComments:Seq=%s\n" % (pep_chr,pep_chr_end,pep_chr_end,peptide))

    #handle exceptions
    if pep_chr_start>pep_chr_end:
        non_mapped_pep+=1
        print(("mapping error",peptide,protein_id, "skip this peptide"))
        continue;
    if pep_chr_start<=0:
        non_mapped_pep+=1
        print(("mapping error",peptide,protein_id,"skip this peptide"))
        continue;

    #print pep_chr_start,pep_chr_end
    #print pep_start_exon,pep_end_exon
    if "chr" not in pep_chr:
        pep_chr="chr"+pep_chr.replace("MT","M") # replace "MT" with "M"

    if pep_start_exon==pep_end_exon: #if peptide map to one exon
        gff_format_line1=[pep_chr,"MS","mRNA",pep_chr_start,pep_chr_end,".",pep_strand,".","ID="+peptide]
        gff_format_line2=[pep_chr,"MS","CDS",pep_chr_start,pep_chr_end,".",pep_strand,"0","Parent="+peptide]
        gff_output.write("\t".join(map(str,gff_format_line1))+"\n")
        gff_output.write("\t".join(map(str,gff_format_line2))+"\n")
    elif abs(pep_start_exon-pep_end_exon)==1: #if it is a splice junction peptide spanning over two exons
        if pep_strand=="+":
            gff_format_line1=[pep_chr,"MS","mRNA",pep_chr_start,pep_chr_end,".",pep_strand,".","ID="+peptide]
            gff_format_line2=[pep_chr,"MS","CDS",pep_chr_start,exons[pep_start_exon-1].end,".",pep_strand,"0","Parent="+peptide]
            gff_format_line3=[pep_chr,"MS","CDS",exons[pep_end_exon-1].start,pep_chr_end,".",pep_strand,".","Parent="+peptide]
        else:
            gff_format_line1=[pep_chr,"MS","mRNA",pep_chr_start,pep_chr_end,".",pep_strand,".","ID="+peptide]
            gff_format_line2=[pep_chr,"MS","CDS",pep_chr_start,exons[pep_end_exon-1].end,".",pep_strand,"0","Parent="+peptide]
            gff_format_line3=[pep_chr,"MS","CDS",exons[pep_start_exon-1].start,pep_chr_end,".",pep_strand,".","Parent="+peptide]

        gff_output.write("\t".join(map(str,gff_format_line1))+"\n")
        gff_output.write("\t".join(map(str,gff_format_line2))+"\n")
        gff_output.write("\t".join(map(str,gff_format_line3))+"\n")
    elif abs(pep_start_exon-pep_end_exon)>1: #if peptide span multiple exons,rare case!
        if pep_strand=="+":
            gff_format_line1=[pep_chr,"MS","mRNA",pep_chr_start,pep_chr_end,".",pep_strand,".","ID="+peptide]
            gff_format_line2=[pep_chr,"MS","CDS",pep_chr_start,exons[pep_start_exon-1].end,".",pep_strand,"0","Parent="+peptide]
            gff_format_line3=[pep_chr,"MS","CDS",exons[pep_end_exon-1].start,pep_chr_end,".",pep_strand,".","Parent="+peptide]
        else:
            gff_format_line1=[pep_chr,"MS","mRNA",pep_chr_start,pep_chr_end,".",pep_strand,".","ID="+peptide]
            gff_format_line2=[pep_chr,"MS","CDS",pep_chr_start,exons[pep_end_exon-1].end,".",pep_strand,"0","Parent="+peptide]
            gff_format_line3=[pep_chr,"MS","CDS",exons[pep_start_exon-1].start,pep_chr_end,".",pep_strand,".","Parent="+peptide]

        gff_output.write("\t".join(map(str,gff_format_line1))+"\n")
        gff_output.write("\t".join(map(str,gff_format_line2))+"\n")
        for k in range(min(pep_start_exon,pep_end_exon)+1,max(pep_start_exon,pep_end_exon)):
            gff_format_line=[pep_chr,"MS","CDS",exons[k-1].start,exons[k-1].end,".",pep_strand,".","Parent="+peptide]
            gff_output.write("\t".join(map(str,gff_format_line))+"\n")

        gff_output.write("\t".join(map(str,gff_format_line3))+"\n")

gff_output.close()
tab_output.close()
bed_output.close()
fasta_output.close()

print(("total number of unique peptides",len(novpep_dic)))
print(("total number of unmapped peptides",non_mapped_pep))
