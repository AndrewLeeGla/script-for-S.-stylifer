#!/usr/bin/python3
#-*- coding: UTF-8 -*-

#by ChaoLi, Apr 13 2019, cleeouc@gmail.com
#usage: python3 contigs_bases_before_genes_analysis.py  -i <genome.fasta>  -I <cds.fasta.id>  -o <genome_contigs_of_cds.fasta>  -m <bases_before_genes_of_5.fasta>  -n <bases_before_genes_of_3.fasta>

import sys, getopt, re, math
from Bio import SeqIO

def main(argv):
    global inputfile
    
    global Inputfile
    
    global outputfile
    
    global outputfile_5
    
    global outputfile_3
    
    global xinfile
    
    try:
        opts, args = getopt.getopt(argv, "hi:I:o:x:m:n:", ["infile=", "Infile=", "xinfile=", "outfile=", "outfile_5=", "outfile_3="])
    except getopt.GetoptError:
        print ('usage: python3 extract_genes_from_gff_and_genome_fasta.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -x <xinfile>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>')
        sys,exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('usage: python3 extract_genes_from_gff_and_genome_fasta.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -x <xinfile>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inputfile = arg
        elif opt in ("-I", "--Infile"):
            Inputfile = arg
        elif opt in ("-o", "--outfile"):
            outputfile = arg
        elif opt in ("-m", "--outfile_5"):
            outputfile_5 = arg
        elif opt in ("-n", "--outfile_3"):
            outputfile_3 = arg
        elif opt in ("-x", "--xinfile"):
            xinfile = arg

def ATGC_percentage_count(infile, outfile):
    a = open(infile, "r")
    b = open(outfile, "w")
    position = {}
    shannon = {}
    
    for i in a.readlines():
        if not i.startswith("position"):
            item = i.rstrip().split("\t")
            position.setdefault(int(item[0]),[])
            position[int(item[0])].append(item[2])
    
    for i in range(len(position)):
        shannon.setdefault(i+1, 0.000)
        for key in position[i+1]:
            shannon[i+1] -= float(key) * math.log(float(key),2)
        shannon[i+1] = round(shannon[i+1],3)
        b.write(str(i+1) + "\t" + "H(U)" + "\t" + str(shannon[i+1]) + "\n")
    b.close()
    
if __name__ == "__main__":
    main(sys.argv[1:])
    ATGC_percentage_count(inputfile, outputfile)