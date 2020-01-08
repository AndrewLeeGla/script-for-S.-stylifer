#!/usr/bin/python3
#-*- coding: UTF-8 -*-

#by ChaoLi, Apr 13 2019, cleeouc@gmail.com
#usage: python3 contigs_bases_before_genes_analysis.py  -i <genome.fasta>  -I <cds.fasta.id>  -o <genome_contigs_of_cds.fasta>  -m <bases_before_genes_of_5.fasta>  -n <bases_before_genes_of_3.fasta>

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
            moutputfile_5 = arg
        elif opt in ("-n", "--outfile_3"):
            noutputfile_3 = arg
        elif opt in ("-x", "--xinfile"):
            xinfile = arg


def extract_genes_from_gff_and_genome_fasta(gff, genome, genes_5, genes_3):
    a = open(gff, encoding = 'utf-8', mode = 'r')
    b = SeqIO.parse(genome, "fasta")
    c = open(genes_5, "w")
    d = open(genes_3, "w")
#    e = open(genes2, "w")
#    f = open(genes3, "w")
    genome_dict= {}
    ################################ 做基因组字典 ####################################
    for m in b:
        genome_dict[m.id] = m
    
    ############################ 提取仅包含一个基因的contig ##########################
    contigs = []
    contigs_gff = []
    for gff_information in a.readlines():
        if not gff_information.startswith("#"):
            items = gff_information.split("\t")
            if items[2] == "transcript":
                contigs_gff.append(items)
                contigs.append(items[0])
    for contig_id in contigs:
        if contigs.count(contig_id) == 1:
        ##############################################################################
            for gff_information in contigs_gff:
                if gff_information[0] == contig_id:
                    if gff_information[6] == "+":
                        c.write(">" + gff_information[0] + "_5_UTR\n" + str(genome_dict[contig_id].seq[: int(items[3])]) + "\n")
                        d.write(">" + gff_information[0] + "_3_UTR\n" + str(genome_dict[contig_id].seq[int(items[4]):]) + "\n")
                    elif gff_information[6] == "-":
                        sequence_before = genome_dict[contig_id].seq[int(gff_information[4]):]
                        c.write(">" + gff_information[0] + "_5_UTR\n" + str(sequence_before.reverse_complement()) + "\n")
                        sequence_behind = genome_dict[contig_id].seq[: int(gff_information[3])]
                        d.write(">" + gff_information[0] + "_3_UTR\n" + str(sequence_behind.reverse_complement()) + "\n")
    a.close()
    b.close()
    c.close()
    d.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    extract_genes_from_gff_and_genome_fasta(inputfile, Inputfile, outputfile, moutputfile_5)