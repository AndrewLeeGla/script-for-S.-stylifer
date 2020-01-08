#!/usr/bin/python3
#-*- coding: UTF-8 -*-

#by ChaoLi, Apr 13 2019, cleeouc@gmail.com
#usage: python3 contigs_bases_before_genes_analysis.py  -i <genome.fasta>  -I <cds.fasta.id>  -o <genome_contigs_of_cds.fasta>  -m <bases_before_genes_of_5.fasta>  -n <bases_before_genes_of_3.fasta>

import sys, getopt, re
from Bio import SeqIO
#from collections import defaultdict

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


def extract_genes_from_gff_and_genome_fasta(blastx, cds, statistic, summary):
    a = open(blastx, encoding = 'utf-8', mode = 'r')
    b = SeqIO.parse(cds, "fasta")
    c = open(statistic, "w")
    d = open(summary, "w")
#    e = open(genes2, "w")
#    f = open(genes3, "w")
    blastx_1 = {}
    cds_fasta = {}
    c.write("stop codon" + "\t" + "amino" + "\n")
    ################################ 做blastx字典 ####################################
    for m in a.readlines():
        if m.startswith("Query="):
            Query_id = m.split()[1]
#            print(Query_id)
            blastx_1[Query_id] = []
            ########################## 解决 **** No hits found **** 的问题 ***********
        elif m.startswith("Query") or m.startswith("Sbjct"):
            blastx_1[Query_id].append(m)
    
    ################################ 做blastx字典 ####################################
    for seq in b:
#        print(seq.id)
        cds_fasta[seq.id] = seq
    
    ############################ 提取stop codon ##########################
    amino_acids = {}
    ############################ 把ami 传入 amino_acids.setdefault(str(stop_codon), ami) 最后出来的字典，虽然键不同，但是每个键的值均为一个 列表 ami，每次根据键，对amino_acids[key]的修改，都是修改的 ami，所以每个键的值都是相同的！################################################
    ami = []
#    TAA = []
#    TAG = []
#    TGA = []
#    contigs_gff = []
    
    for key in blastx_1:
        if blastx_1[key]:
            for line in blastx_1[key]:
#            print("#####################################" + "\n" + line[2])
                lines = line.split()
                if re.search('\*',lines[2]):
                    f = re.finditer('\*',lines[2])
                    for i in f:
    #                    print(i)
                        pos = int(i.span()[0])
                        ####################### 判断*前后的氨基酸种类一致 #####################################
                        if (int(pos-1) >= 0) and (len(lines[2])-1 >= int(pos+1)):
    #                        print(str(lines) + "\t" + str(pos))
    #                        print(i.span()[1])
    #                        print(lines[2] + "##############" + str(len(lines[2])))
                            subject = blastx_1[key][blastx_1[key].index(line)+1].split()
                            if (lines[2][pos-1] == subject[2][pos-1]) and (lines[2][pos+1] == subject[2][pos+1]):
                                print(lines)
                                print(i.span()[0])
                                print(lines[2])
                                print(subject[2])
                                print(str(cds_fasta[key].seq))
                                if int(lines[1]) < int(lines[3]):
                                    stop_codon = str(cds_fasta[key].seq[int(int(lines[1]) + (pos * 3) - 1):int(int(lines[1]) + (pos*3) + 2)])
    #                                stop_codon = str(cds_fasta[key].seq)[int(lines[1])+pos*3-1:int(lines[1])+ pos*3 + 2]
                                    amino_acids.setdefault(str(stop_codon), [])
                                    amino_acids[str(stop_codon)].append(subject[2][pos])
#                                    amino_acids[str(stop_codon)] = amino_acids[str(stop_codon)]
                                    c.write(stop_codon + "\t" + subject[2][pos] + "\n")
                                elif int(lines[1]) > int(lines[3]):
    #                                amino_acids[stop_codon].append(subject[2][pos])
                                    stop_codon = str(cds_fasta[key].seq[int(int(lines[1]) - (pos * 3) - 2 - 1):int(int(lines[1]) - (pos*3))].reverse_complement())
                                    amino_acids.setdefault(str(stop_codon), [])
                                    amino_acids[str(stop_codon)].append(subject[2][pos])
#                                    amino_acids[str(stop_codon)] = amino_acids[str(stop_codon)]
    #                                stop_codon = str(cds_fasta[key].seq)[int(lines[1])-pos*3-2-1:int(lines[1])-pos*3]
                                    c.write(stop_codon + "\t" + subject[2][pos] + "\n")
    for k in amino_acids:
#        print(k)
#        print(v)
        print(str(k), amino_acids[str(k)])
        d.write("stop codon " + str(k) + " number: " + str(len(amino_acids[k])) + "\n")
#    for key in amino_acids:
    for k in amino_acids:
        don = []
#        print(key)
#        print(amino_acids[key])
        if str(k) in ["TAA", "TAG", "TGA"]:
#            print(amino_acids[key])
            for i in amino_acids[k]:
                if i not in don:
                    don.append(i)
                    d.write("coded by " + str(k) + " the " + i + " number is: " + str(amino_acids[k].count(i)) + "\n")
    a.close()
    b.close()
    c.close()
    d.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    extract_genes_from_gff_and_genome_fasta(inputfile, Inputfile, outputfile, moutputfile_5)