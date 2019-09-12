#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys

syntax = '''
------------------------------------------------------------------------------------
Usage: python interleave_fq2fa.py R1.fastq R2.fastq output
------------------------------------------------------------------------------------
'''

if len(sys.argv) != 4:
    print syntax
    sys.exit()

outname = sys.argv[3]


#--Reading files 
infile1 = open (sys.argv[1], 'r')
infile2 = open (sys.argv[2], 'r')
outfile = open (outname, 'w')
r = 0

while True:
                
    name1 = infile1.readline().rstrip('\n')
    seq1 = infile1.readline().rstrip('\n')
    coment1 = infile1.readline().rstrip('\n')
    qual1 = infile1.readline().rstrip('\n')

    name2 = infile2.readline().rstrip('\n')
    seq2 = infile2.readline().rstrip('\n')
    coment2 = infile2.readline().rstrip('\n')
    qual2 = infile2.readline().rstrip('\n')

    if not name1:
        break

    if name1.split()[0] != name2.split()[0]:
        print 'ERROR: reads are not properly sorted'
        print name1, name2
        break

    outfile.write('>read' + str(r) + '/1' + '\n')
    outfile.write(seq1 + '\n')
    outfile.write('>read' + str(r) + '/2' + '\n')
    outfile.write(seq2 + '\n')
    
    r = r + 1
    
infile1.close()
infile2.close()
outfile.close()

