#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys


##########################################################################################
syntax = '''
------------------------------------------------------------------------------------
Usage: python rename_contigs.py contigs_file.fa 
------------------------------------------------------------------------------------
'''
##########################################################################################

if len(sys.argv) != 2:
        print syntax
        sys.exit()

##########################################################################################

#--Variables
fasta = {}
contigs = []
prefix = sys.argv[1].split('.')[0]

#--We screen the file that contains the contigs in fasta format
infile = open (sys.argv[1], 'r')
seq = ""
for line in infile:
        line = line.rstrip('\n')

        if line[0] == '>':
                if seq:
                        fasta[name] = seq
                        seq = ""
                name = line[1:]
                contigs.append(name)
        else:
                seq = seq + line

#--The last contig
fasta[name] = seq
infile.close()

#--We write the sequences

outfile = open (prefix + '_renamed.fa', 'w')
i = 0
for name in contigs:

        new_name = '>' + prefix + '_' + str(i) + ' ' + name
        seq = fasta[name]

        outfile.write(new_name + '\n')
        outfile.write(seq + '\n')
        i = i + 1

outfile.close()

