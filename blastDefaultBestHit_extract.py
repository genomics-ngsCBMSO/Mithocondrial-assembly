#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import os
import argparse
from collections import defaultdict

##########################################################################################

if __name__ == '__main__':
    
    ##########################################################################################
    #--Argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest = 'blast', type = str, required = True, help = 'blast_file outfmt default')
    args = parser.parse_args()
    ##########################################################################################
    
    blast_file = open(args.blast, 'r')

    contigs = {}
    t = 0
    for line in blast_file:

        line = line.rstrip('\n')

        if not line: # to avoid empty lines
            continue

        if 'Query=' in line: # contig name
            best_hit = ''
            contig = line.split()[1]
            t = 0  
            continue

        if "***** No hits found *****" in line: # no hit was found for this query
            continue

        if "Sequences producing significant alignments:" in line:
            t = 1
            continue

        if t == 1:
            hit = line.split()
            eval = float(hit[-1])
            score = float(hit[-2])
            
            if eval <= 1e-3 or score >= 35:
                t = 2
            else:
                t = 0
                continue

            continue
        
        if t == 2:
            
            if line[0] == '>':
                best_hit = line[1:]
                t = 3
                continue
                
        if t == 3:
            best_hit = best_hit + ' ' + line
            if 'Length=' in line:
                contigs[contig] = [best_hit, str(eval), str(score)]
                t = 0
            continue

    output = open(args.blast.split('.')[0] + '_table.txt', 'w')
    for contig, hit in contigs.iteritems():
        output.write(contig + '\t' + '\t'.join(hit) + '\n')

    output.close()
