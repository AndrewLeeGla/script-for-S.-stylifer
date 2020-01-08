#!/usr/bin/python3
#-*- coding: UTF-8 -*-

#by ChaoLi, May 23 2019, cleeouc@gmail.com
#usage: python3 Stylifer_telomere_analysis.py  -i <genome.fasta>  -I <cds.fasta.id>  -o <genome_contigs_of_cds.fasta>  -m <bases_before_genes_of_5.fasta>  -n <bases_before_genes_of_3.fasta>


import re
import sys, getopt
from typing import TextIO
from Bio import SeqIO

def main(argv):
    global inputfile
    
    global Inputfile
    
    global outputfile
    
    global outputfile_5
    
    global outputfile_3
    
    global xinfile
    
    global poutfile
    
    global qoutfile
    
    global routfile
    
    global soutfile
    
    global toutfile
    
    
    try:
        opts, args = getopt.getopt(argv, "p:q:r:s:t:h:i:I:o:x:m:n:x:", ["infile=", "Infile=", "xinfile=", "outfile=", "outfile_5=", "outfile_3="])
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
        elif opt in ("-p"):
            poutfile = arg
        elif opt in ("-q"):
            qoutfile = arg
        elif opt in ("-r"):
            routfile = arg
        elif opt in ("-s"):
            soutfile = arg
        elif opt in ("-t"):
            toutfile = arg

#out1 = open(r"/home/zhengwb/Strombidium_stylifer/02_telomere_analysis/telo2.fasta", "w")  # double telomere
#out2 = open(r"/home/zhengwb/Strombidium_stylifer/02_telomere_analysis/telo1.fasta", "w")  # single telomere
#out3 = open(r"/home/zhengwb/Strombidium_stylifer/02_telomere_analysis/telo0.fasta", "w")  # zero telomere
#out4 = open(r"/home/zhengwb/Strombidium_stylifer/02_telomere_analysis/nothing.fa", "w")  # waste
#out5 = open(r"/home/zhengwb/Strombidium_stylifer/02_telomere_analysis/no_redundance.fasta", "w")  # combine
#out6 = open(r"/home/zhengwb/Strombidium_stylifer/02_telomere_analysis/0telo_highGC.fa", "w")  # highGC
def telomere_analysis(input, out1n, out2n, out3n, out4n, out5n, out6n):
    out1 = open(out1n, "w")
    out2 = open(out2n, "w")
    out3 = open(out3n, "w")
    out4 = open(out4n, "w")
    out5 = open(out5n, "w")
    out6 = open(out6n, "w")

    for rec in SeqIO.parse(input, "fasta"):
        seq_this = str(rec.seq)
        OK5 = 0
        OK3 = 0
        length_seq = len(seq_this)
        GC = (seq_this.count("G") + seq_this.count("C")) / length_seq
        mobj1 = re.search("CCCCAAAACCCC", seq_this)
        mobj2 = re.search("GGGGTTTTGGGG", seq_this)
        tag = rec.id
        if mobj1:
            OK5 = 1
        if mobj2:
            OK3 = 1
        if OK5 == 1 and OK3 == 1 and length_seq >= 200:
            for m in re.finditer("CCCCAAAACCCC", seq_this):
                pos_1 = m.start()
            for n in re.finditer("GGGGTTTTGGGG", seq_this):
                pos_2 = n.start()
            distance = abs(pos_2 - pos_1)
            if distance > 200:
                print(">%s\n%s" % (tag, seq_this), file = out1)
                print(">%s\n%s" % (tag, seq_this), file = out5)
            else:
                print(">%s\n%s" % (tag, seq_this), file = out4)
        if OK5 == 1 and OK3 != 1 and length_seq >= 200:
            print(">%s\n%s" % (tag, seq_this), file = out2)
            print(">%s\n%s" % (tag, seq_this), file=out5)
        if OK5 != 1 and OK3 == 1 and length_seq >= 200:
            print(">%s\n%s" % (tag, seq_this), file = out2)
            print(">%s\n%s" % (tag, seq_this), file = out5)
        if OK3 != 1 and OK5 != 1:
            if length_seq >= 1000:
                if GC < 0.5:
                    print(">%s\n%s" % (tag, seq_this), file = out3)
                    print(">%s\n%s" % (tag, seq_this), file = out5)
                else:
                    print(">%s\n%s" % (tag, seq_this), file = out6)
            else:
                print(">%s\n%s" % (tag, seq_this), file = out4)

if __name__ == "__main__":
    main(sys.argv[1:])
    telomere_analysis(inputfile, outputfile, poutfile, qoutfile, routfile, soutfile, toutfile)