#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import os
import argparse
from collections import defaultdict


##########################################################################################

def checkfile (file): 

    """ File checking """

    if not os.path.exists(file):
        print "-" * 90
        print "ERROR: file '" + file + "' not found"
        print "-" * 90
        sys.exit()

def revcomp(seq):

    """ This function returns the reverse-complement of a sequence """

    seq = seq[::-1].upper()
    basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N' : 'N', 'R':'Y', 'Y':'R', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'N':'N'}
    seq = list(seq) 
    seq = [basecomplement[base] for base in seq] 
    return ''.join(seq) 

def fasta2dict(file):
    
    """ This function reads a fasta file and stores the sequences in a dictionary """
    
    #--Variables
    fasta = {}
    seq = "" 

    infile = open (file, 'r')
    
    for line in infile:
        line = line.rstrip('\n')

        if line[0] == '>':
        
            if seq:
            
                fasta[name] = seq
                seq = ""
                
            name = line.split()[0][1:]
            
        else:
        
            seq = seq + line

    #--Last sequence
    fasta[name] = seq
    infile.close()
    
    #--Return dictionary object
    return fasta

##########################################################################################

if __name__ == '__main__':
    
    ##########################################################################################
    #--Argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest = 'fasta_file', type = str, required = True, help = 'multifasta or single fasta file')
    parser.add_argument('-l', dest = 'min_prot_len',  type = int, default = 80, help = 'Minimun protein length')
    args = parser.parse_args()
    ##########################################################################################
    
    #--Variables
    min_prot_len = args.min_prot_len
    
    #--Checking input file
    checkfile(args.fasta_file)

    #--Reading input file
    fasta = fasta2dict(args.fasta_file)

    #--Generating genetic code translation table (Protozoan Mitochondrial Code)

    # source = http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG5
    # The Standard Code
    #   AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    # Starts = ---M---------------M---------------M----------------------------
    # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    #   AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    # Starts = --MM---------------M------------MMMM---------------M------------

    aas  =   "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    start_codons = ['ATG', 'ATT', 'ATA'] # Leishmania start codons
    Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
    
    code = {}
    for i, b1 in enumerate(Base1):
        aa = aas[i]
        b2 = Base2[i]
        b3 = Base3[i]
        code[b1 + b2 + b3] = aa

    #--Analyzing seqs
    results = {}

    for name, seq in fasta.iteritems():

        seqlen = len(seq)

        for strand in ['+', '-']: #--Each strand 

            if strand == '+':
                temp_seq = seq.upper()
            else:
                temp_seq = revcomp(seq.upper()) 

            for frame in range(3):  #--Each frame
            
                codons = [temp_seq[i:i+3] for i in xrange(frame,len(temp_seq), 3)] #--Split sequence in codons
                prot = [code.get(codon, 'X') for codon in codons]  #--Translate complete sequece
                stops = [i for i, p in enumerate(prot) if p == '*'] #--Finding stops in the protein sequence

                #--Analysing first peptide. This peptide could not have a start codon, but it could be upstream, although we do not know the upstream sequence.
                stop = stops[0]
                peptide = prot[0: stop]

                if len(peptide) >= min_prot_len: #--Length pre-filter

                        starts = [i for i, p in enumerate(peptide) if codons[i + start] in start_codons]

                        if starts:

                            if len(peptide[starts[0]:]) >= min_prot_len: #--Length filter using the appropiate start codon

                                if strand == '+':
                                    x = (starts[0] + start) * 3 + frame
                                    y = (stop + 1) * 3 + frame
                                    results[(x, y, name)] = [name, x, y, strand, ''.join(peptide[starts[0]:]), 'complete', seq]


                                else:
                                    x = (seqlen - ((stop + 1) * 3 + frame)) + 1
                                    y = seqlen - ((starts[0] + start) * 3 + frame)
                                    results[(x, y, name)] = [name, x, y, strand, ''.join(peptide[starts[0]:]), 'complete', seq]

                        else:

                            if strand == '+':
                                x = frame + 1
                                y = (stop + 1) * 3 + frame
                                results[(x, y, name)] = [name, x, y, strand, ''.join(peptide), 'partial', seq]

                            else:
                                x = (seqlen - ((stop + 1) * 3 + frame))
                                y = seqlen - (frame + 1)
                                results[(x, y, name)] = [name, x, y, strand, ''.join(peptide), 'partial', seq]

                #--Analysing inner peptides
                start = stops[0] 
                for stop in stops[1:]:

                    peptide = prot[start: stop]

                    if len(peptide) >= min_prot_len: #--Length pre-filter

                        starts = [i for i, p in enumerate(peptide) if codons[i + start] in start_codons]
            
                        if starts:

                            if len(peptide[starts[0]:]) >= min_prot_len: #--Length filter using the appropiate start codon

                                if strand == '+':
                                    x = (starts[0] + start) * 3 + frame
                                    y = (stop + 1) * 3 + frame
                                    results[(x, y, name)] = [name, x, y, strand, ''.join(peptide[starts[0]:]), 'complete', seq]

                                else:
                                    x = (seqlen - ((stop + 1) * 3 + frame)) + 1
                                    y = seqlen - ((starts[0] + start) * 3 + frame)
                                    results[(x, y, name)] = [name, x, y, strand, ''.join(peptide[starts[0]:]), 'complete', seq]
                    
                    start = stop + 1

                #--Analysing last peptide. This peptide do not have stop, but could be downstream, although we do not know the downstream sequence.
                start = stops[-1]
                stop = len(codons) + 1
                peptide = prot[start: stop]

                if len(peptide) >= min_prot_len: #--Length pre-filter

                        starts = [i for i, p in enumerate(peptide) if codons[i + start] in start_codons]
            
                        if starts:

                            if len(peptide[starts[0]:]) >= min_prot_len: #--Length filter using the appropiate start codon

                                if strand == '+':
                                    x = (starts[0] + start) * 3 + frame
                                    y = (stop + 1) * 3 + frame
                                    results[(x, y, name)] = [name, x, y, strand, ''.join(peptide[starts[0]:]), 'partial', seq]

                                else:
                                    x = (seqlen - ((stop + 1) * 3 + frame)) + 1
                                    y = seqlen - ((starts[0] + start) * 3 + frame)
                                    results[(x, y, name)] = [name, x, y, strand, ''.join(peptide[starts[0]:]), 'partial', seq]

    #--Writing results
    prefix = args.fasta_file.split('.')[0]
    prots_file = open(prefix + '_prots.fa', 'w')
    orfs_file = open(prefix + '_orfs.fa', 'w')
    gtf_file = open(prefix + '_orfs.gtf', 'w')
    
    c = 1
    for orf in sorted(results.keys()):
        name, x, y, strand, peptide, type, seq = results[orf]
        if strand == '+':
            locus = name + ':' + str(x + 1) + '-' + str(y) + ':' + strand
        else:
            locus = name + ':' + str(x) + '-' + str(y) + ':' + strand
        
        if type == 'partial':
            new_name = name + '_partialORF_' + str(c)
        else:
            new_name = name + '_ORF_' + str(c)
        
        
        #--Protein sequences
        prots_file.write('>' + new_name + ' ' + locus + '\n')
        prots_file.write(peptide + '\n')
        
        #--ORF sequences
        orfs_file.write('>' + new_name + ' ' + locus + '\n')
        if strand == '+':
            orfs_file.write(seq[x:y] + '\n')
        else:
            orfs_file.write(revcomp(seq[x-1:y]) + '\n')

        #--GTF file
        if strand == '+':
            col = [name, 'CBMSO', 'gene', str(x + 1), str(y), '.', strand, '.']
        else:
            col = [name, 'CBMSO', 'gene', str(x), str(y), '.', strand, '.']

        att = ['gene_id "',new_name, '"; name "', new_name, '"; coordinates "', locus, '";']
        gtf_file.write('\t'.join(col) + '\t' + ''.join(att) + '\n')     

        c = c + 1


    #--Closing files
    gtf_file.close()
    prots_file.close()
    orfs_file.close()
