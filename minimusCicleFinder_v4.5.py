#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from __future__ import division
import sys
import subprocess
import os
import argparse

def fasta2dict(file):
    
    """ This fuction screens a fasta file and stock the information within a dictionary """
    
    #--Variables
    fasta = {}
    seq = "" 

    #-- We screen the file that contains the fasta sequences 
    infile = open (file, 'r')
    
    for line in infile:
        line = line.rstrip('\n')

        if line[0] == '>':
        
            if seq:
            
                fasta[name] = seq
                seq = ""
                
            name = line.split()[0]
            
        else:
        
            seq = seq + line

    #--The last sequence
    fasta[name] = seq
    infile.close()
    
    #--We return the dictionary
    return fasta
    
def fastaseq(file):
    
    fasta_seqs = []

    seq = ''
    name = ''

    #--We screen the file that contains the fasta sequences 
    infile = open (file, 'r')

    for line in infile:
        line = line.rstrip('\n')

        if line[0] == '>':
            continue
        else:
            seq = seq + line

    infile.close()
    
    #--We return the dictionary
    return seq
    
##########################################################################################
#--Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest = 'input', type = str, required = True, help = 'Input file in fasta format')
parser.add_argument('-o', dest = 'overlap', type = int, default = 40, help = 'Minimum overlap, default 40')
parser.add_argument('-r', dest = 'refcount', type = int, default = 1, help = 'Number of sequences is the first set, (Def 1)')
parser.add_argument('-c', dest = 'conserr', type = float, default = 0.06, help = 'Maximum consensus error (0..1) (Def 0.06)')
parser.add_argument('-m', dest = 'minid', type = float, default = 94.0, help = 'Minimum overlap percentage identity for align. (Def 94)')
parser.add_argument('-t', dest = 'maxtrim', type = int, default = 20, help = 'Maximum sequence trimming length (Def 20bp)')
parser.add_argument('-k', dest = 'kmer', type = int, required = True, help = 'max kmer used in the assembly')
args = parser.parse_args()

##########################################################################################
#--Parameters and varibles  
input = args.input
if not os.path.exists(input):
    print 'ERROR: file ' + file + ' does not exist'
    sys.exit()

kmer = args.kmer
prefix = input.split('.')[0]

##########################################################################################
#--Reading file
fasta = fasta2dict(input)
names = sorted(fasta.keys())
c = 0

#--We delete the former temporary files
subprocess.call('rm -fr lr.*', shell = True)
subprocess.call('rm -fr right.fa left.fa', shell = True) 


#--We create the output file
output = open(prefix + '_circles.txt', 'w')

for name in names:
    seq = fasta[name]
    clen = len(seq) 

    #--We create two files so that they contain the analysed sequences
    outR = open('right.fa', 'w')
    outL = open('left.fa', 'w')
    seql = seq[:len(seq)//2]
    seqr = seq[len(seq)//2:]
    outL.write(name + '_L' + '\n' + seql + '\n')
    outR.write(name + '_R' + '\n' + seqr + '\n')
    outR.close()
    outL.close()
    
    #--We run minimus2 (it has to be in the same folder as the script so that it works properly)
    log = open('temp.log', 'w')
    subprocess.call('cat left.fa right.fa > lr.fa', shell = True, stdout=log, stderr=log) 
    subprocess.call('toAmos -s lr.fa -o lr.afg', shell = True, stdout=log, stderr=log) 
    subprocess.call('/PATHtoWD/minimus2 lr -D REFCOUNT=' + str(args.refcount) + ' -D OVERLAP=' + str(args.overlap) + ' -D CONSERR=' + str(args.conserr) + ' -D MINID=' + str(args.minid) + ' -D MAXTRIM=' + str(args.maxtrim), shell = True, stdout=log, stderr=log)
    log.close()
    
    #--We corroborate that the output file has been generated
    if os.path.exists('lr.fasta'):
    
        #-- Number of sequences in output fasta file control
        out, err = subprocess.Popen(["grep", "-c", ">", 'lr.fasta'], stdout=subprocess.PIPE).communicate()
        
        if int(out) == 0: #--The output file was created but not a proper overlap accordingly to the defined parameters
            output.write(name + ' length=' + str(len(seq)) + '\n' + seq + '\n')
            
        elif int(out) == 1: #--Overlap was found
            circle = fastaseq('lr.fasta')
            lendiff = len(seq) - len(circle)
            if lendiff > kmer * 1.5:
                output.write(name + ' length=' + str(len(seq)) + '\n' + seq + '\n')
            else: #--Discard longer overlaps than kmer * 1.5
                output.write(name + '_circle length=' + str(len(seq)) + ' lenDiff=' + str(lendiff) + '\n' + circle + '\n')
                
            c = c + 1
        else: 
            print 'ERROR: more than 1 sequence found in lr.fasta'
            subprocess.call('mv lr.fa ' + name + '.fa', shell = True)
            subprocess.call('mv lr.coords ' + name + '.coords', shell = True)
            output.write(name + ' length=' + str(len(seq)) + '\n' + seq + '\n')
        
    else: #--No overlap was found
        output.write(name + ' length=' + str(len(seq)) + '\n' + seq + '\n')
            

    #--We delete the former temporary files
    subprocess.call('rm -fr lr.*', shell = True)
    subprocess.call('rm -fr right.fa left.fa', shell = True) 

output.close()

#--Summary
print 'There are ' + str(c) + ' circular contigs contigs'
