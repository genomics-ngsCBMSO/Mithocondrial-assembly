#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys

syntax = '''
------------------------------------------------------------------------------------
Usage: python rePaired.py R1.fastq R2.fastq

Description: this script reformat paired-end reads file that have been filtered, 
returning 3 fastq files, two for paired reads and one with orfan reads.
------------------------------------------------------------------------------------
'''

if len(sys.argv) != 3:
    print syntax
    sys.exit()


#--Reading file 1
infile = open (sys.argv[1], 'r')
r1 = {}

while True:
                
    name = infile.readline().rstrip('\n')
    seq = infile.readline().rstrip('\n')
    coment = infile.readline().rstrip('\n')
    qual = infile.readline().rstrip('\n')

    if not name:
        break

    r1[name.split()[0]] = [name, seq, coment, qual]

infile.close()

#--Reading file 2
infile = open (sys.argv[2], 'r')
pair1 = open (sys.argv[1].split('.')[0] + '_paired.fastq', 'w')
pair2 = open (sys.argv[2].split('.')[0] + '_paired.fastq', 'w')
orfan = open (sys.argv[1].split('_')[0] + '_orfan.fastq', 'w')

while True:
    name = infile.readline().rstrip('\n')
    seq = infile.readline().rstrip('\n')
    coment = infile.readline().rstrip('\n')
    qual = infile.readline().rstrip('\n')
    record = [name, seq, coment, qual]

    if not name:
        break

    try:
        pair1.write('\n'.join(r1[name.split()[0]]) + '\n')
        pair2.write('\n'.join(record) + '\n')  
        del r1[name.split()[0]]

    except KeyError:
        orfan.write('\n'.join(record) + '\n')
        
#--Writting orfan reads
for read, record in r1.iteritems():
    orfan.write('\n'.join(record) + '\n')


infile.close()
pair1.close()
pair2.close()
orfan.close()
