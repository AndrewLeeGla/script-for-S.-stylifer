#!/usr/bin/python3
#-*- coding: UTF-8 -*-

#by ChaoLi, Apr 13 2019, cleeouc@gmail.com
#usage: python3 find_AT_most_length_and_pos_in_first_200bp_of_5_and_3_region.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>

import sys, getopt, re
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
    
    try:
        opts, args = getopt.getopt(argv, "hi:I:o:x:m:n:p:q:r:s:t:", ["infile=", "Infile=", "xinfile=", "outfile=", "outfile_5=", "outfile_3="])
    except getopt.GetoptError:
        print ('usage: python3 find_AT_most_length_and_pos_in_first_200bp_of_5_and_3_region.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -x <xinfile>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>')
        sys,exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('usage: python3 find_AT_most_length_and_pos_in_first_200bp_of_5_and_3_region.py  -i <genome.fasta>  -I <cds.fasta>  -o <genome_contigs_of_cds.fasta>  -x <xinfile>  -m <first_500bp_of_5.fasta>  -n <first_500bp_of_3.fasta>')
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
            xoutfile = arg
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


def find_AT_most_length_and_pos_in_first_200bp_of_5_and_3_region(fasta_5, fasta_3, output_5UTR, poutfile_3UTR):
    a = SeqIO.parse(fasta_5, "fasta")
    b = SeqIO.parse(fasta_3, "fasta")
    c = open(output_5UTR, "w")
    d = open(poutfile_3UTR, "w")
    
    c.write("region" + "\t" + "AT_size" + "\t" + "position" + "\n")
    d.write("region" + "\t" + "AT_size" + "\t" + "position" + "\n")
    
#    AT_frame = {} #键为每一条contig中最长ATframe的长度；值为该长度ATframe的位置
    
    for i in a:
        pos = 0
        AT_id = 1 #longest_AT的键
#        length = 0 #longest_AT的值
        length_AT = {} #键为每一条contig中每一个ATframe的编号；值为该编号ATframe的长度
        pos_AT = {} #键为每一条contig中每一个ATframe的编号；值为该编号ATframe的位置
        
        for thisbase in str(i.seq):
            if pos == 0:
                former = thisbase
            if pos != 0:
                if (thisbase == "A" or thisbase == "T") and (former == "A" or former == "T"):
                    length_AT[AT_id] = length_AT.setdefault(AT_id, 1) + 1
                    pos_AT[AT_id] = len(i) - pos + 1
                    
                elif (thisbase == "A" or thisbase == "T") and (former == "C" or former == "G"):
                    AT_id = AT_id + 1
                    
            former = thisbase
            pos = pos + 1
#        sorted_length_AT = sorted(length_AT.items(), key = lambda item:item[1], reverse = True)
        print(length_AT)
        print("-------------------------------------")
#        AT_frame[sorted_length_AT[0][0]] = AT_frame.setdefault(sorted_length_AT[0][0], 0) + 1
        for AT in length_AT.items():
        
            if int(AT[1]) >= 5:
                c.write("5 primer region" + "\t" + str(AT[1]) + "\t" + str(pos_AT[AT[0]]) + "\n")

    
    print("----------------------------------------------------------------------------------")
    
    for i in b:
        pos = 0
        AT_id = 1 #longest_AT的键
#        length = 0 #longest_AT的值
        length_AT = {} #键为每一条contig中每一个ATframe的编号；值为该编号ATframe的长度
        pos_AT = {} #键为每一条contig中每一个ATframe的编号；值为该编号ATframe的位置
        
        
        for thisbase in str(i.seq):
            if pos == 0:
                former = thisbase
            if pos != 0:
                if (thisbase == "A" or thisbase == "T") and (former == "A" or former == "T"):
                    length_AT[AT_id] = length_AT.setdefault(AT_id, 1) + 1
                    pos_AT[AT_id] = pos - int(length_AT[AT_id]) + 1
                    
                elif (thisbase == "A" or thisbase == "T") and (former == "C" or former == "G"):
                    AT_id = AT_id + 1
                    
            former = thisbase
            pos = pos + 1
#        sorted_length_AT = sorted(length_AT.items(), key = lambda item:item[1], reverse = True)
        print(length_AT)
        print("-------------------------------------")
#        AT_frame[sorted_length_AT[0][0]] = AT_frame.setdefault(sorted_length_AT[0][0], 0) + 1
        for AT in length_AT.items():
        
            if int(AT[1]) >= 5:
                d.write("3 primer region" + "\t" + str(AT[1]) + "\t" + str(pos_AT[AT[0]]) + "\n")
    c.close()
    d.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    find_AT_most_length_and_pos_in_first_200bp_of_5_and_3_region(inputfile, Inputfile, outputfile, poutfile)