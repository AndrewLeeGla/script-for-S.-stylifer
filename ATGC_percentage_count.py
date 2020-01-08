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
#    a = open(infile, "r")
    b = open(outfile, "w")
#    m = []
    n = {}
    number = {}
    bins = 250
    b.write("position" + "\t" + "type" + "\t" + "percentage" + "\n")
    for i in SeqIO.parse(infile, "fasta"):
        n[i.id] = {"A": {}, "T": {}, "G": {}, "C": {}}
        
        binsize = round((len(i)/bins),5)
        
        for m in range(1,bins + 1):
            n[i.id]["A"][m] = float(i.seq[int((m-1)*binsize):int(m*binsize)].count("A")/binsize)
            n[i.id]["T"][m] = float(i.seq[int((m-1)*binsize):int(m*binsize)].count("T")/binsize)
            n[i.id]["G"][m] = float(i.seq[int((m-1)*binsize):int(m*binsize)].count("G")/binsize)
            n[i.id]["C"][m] = float(i.seq[int((m-1)*binsize):int(m*binsize)].count("C")/binsize)
#        key = (int(float(i.id.split('_')[5])/barsize)+1)*barsize
#        n[key] = n.setdefault(key, 0) + len(i)
#        number[key] = number.setdefault(key, 0) + 1
#        
#    keys = list(n.keys())
#    keys.sort()
#    print(keys)

    type = ["A", "T", "G", "C"]
    ids = list(n.keys())
    for base in type:
        for w in range(1,bins + 1):
            sum = 0
#            float(sum)
            for i in ids:
#                print(i)
                sum = sum + n[i][base][w]
            average = round(sum/len(ids),2)
            b.write(str(w) + "\t" + base + "\t" + str(average) + "\n")
    
    
#    SeqIO.write(m, outfile, "fasta")
#            c = open("temp", "r")
    b.close()

if __name__ == "__main__":
    main(sys.argv[1:])
    ATGC_percentage_count(inputfile, outputfile)