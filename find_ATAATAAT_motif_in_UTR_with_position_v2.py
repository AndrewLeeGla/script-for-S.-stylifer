#!/usr/bin/python3
#-*- coding: UTF-8 -*-

#by ChaoLi, Apr 13 2019, cleeouc@gmail.com
#usage: python3 rm_telomeres_in_one_gene_one_contig.py  -i <genome.fasta>  -I <cds.fasta.id>  -o <genome_contigs_of_cds.fasta>  -m <bases_before_genes_of_5.fasta>  -n <bases_before_genes_of_3.fasta>

import sys, getopt, re
from Bio import SeqIO


def main(argv):
    global inputfile
    
    global Inputfile
    
    global outputfile
    
    global moutputfile_5
    
    global noutputfile_3
    
    global xinfile
    
    try:
        opts, args = getopt.getopt(argv, "hi:I:o:x:m:n:", ["infile=", "Infile=", "xinfile=", "outfile=", "outfile_5=", "outfile_3="])
    except getopt.GetoptError:
        print ('usage: python3 rm_telomeres_in_one_gene_one_contig.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -x <xinfile>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>')
        sys,exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('usage: python3 rm_telomeres_in_one_gene_one_contig.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -x <xinfile>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inputfile = arg
        elif opt in ("-I", "--Infile"):
            Inputfile = arg
        elif opt in ("-o", "--outfile"):
            outputfile = arg
        elif opt in ("-m", "--outfile_5"):
            moutputfile_5 = arg
        elif opt in ("-n", "--outfile_3"):
            noutputfile_3 = arg
        elif opt in ("-x", "--xinfile"):
            xinfile = arg

def rm_telomeres_in_one_gene_one_contig(file_5, file_3, outfile5, AT_2_5):
#, file_3, outfile5, outfile3, AT_2_5, AT_2_3
    a = SeqIO.parse(file_5, "fasta")
    b = SeqIO.parse(file_3, "fasta")
    c = open(outfile5, "w")
#    d = open(outfile3, "w")
    e = open(AT_2_5, "w")
#    f = open(AT_2_3, "w")
    
    AT_1 = re.compile(r'AT[AT]AT[AT]A[AT]')
#    AT_2 = re.compile(r'T.AT.ATTATT')
    
    c.write("region" + "\t" + "pattern" + "\t" + "position" + "\n")
#    d.write("region" + "\t" + "pattern" + "\t" + "position" + "\n")
    e.write("region" + "\t" + "pattern" + "\t" + "position" + "\n")
#    f.write("region" + "\t" + "pattern" + "\t" + "position" + "\n")
    
    for i in a:
        if AT_1.findall(str(i.seq)):
            for AT in AT_1.finditer(str(i.seq)):
                c.write("5_UTR" + "\t" + str(AT.group()) + "\t" + str(len(i) - (AT.start() + 1) + 1) + "\n")

#####
    for i in b:
        if AT_1.findall(str(i.seq)):
            for AT in AT_1.finditer(str(i.seq)):
                e.write("3_UTR" + "\t" + str(AT.group()) + "\t" + str(AT.start() + 1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
    rm_telomeres_in_one_gene_one_contig(inputfile, Inputfile, moutputfile_5, noutputfile_3)
    #, Inputfile, moutputfile_5, outputfile, noutputfile_3, xinfile