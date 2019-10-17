# Mitochondrial Assembly


## INTRODUCTION

This pipeline can be used in order to assembly the maxicircle and the minicircles
of kinetoplastids from DNA-sequencing data. It has been successfully tested on different 
*Leishmania* species: *L. major*, *L. infantum*, *L. braziliensis*, *L. adleri* and *L. guyanensis*.

## REQUIREMENTS

This pipeline requires the following dependences: 

* [python2](https://www.python.org/downloads/)
* [perl](https://www.perl.org/)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://samtools.sourceforge.net/)
* [Bedtools](https://bedtools.readthedocs.io/en/latest/)
* [IDBA_UD](https://github.com/loneknightpy/idba)
* [toAmos](https://sourceforge.net/projects/amos/)  (which includes Minimus2)
* [ORF-Finder](https://github.com/averissimo/orf_finder)
* [Blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [prinseq](http://prinseq.sourceforge.net/)


All these softwares must be in the local PATH except Minimus2 (see minimusCicleFinder_v4.5.py for further information).


## INSTALLATION

Make sure you have installed all the dependencies before try to run the pipeline!

The downloaded folder contains  all the *in-house* additional scripts needed to do the complete pipeline:

* rePaired.py
* interleave_fq2fa.py
* rename_contigs.py
* blastDefaultBestHit_extract.py
* minimusCicleFinder_v4.5.py


## USAGE


### 1. Alignment with Bowtie2

Each pair is aligned separately, using the following parameters:

```
bowtie2 --local -p 8 -x <index> -U <forward read.fastq> -S <alignment_forward.sam>
bowtie2 --local -p 8 -x <index> -U <reverse_read.fastq> -S <alignment_reverse.sam>

```

This step will generate 2 files:

- alignment_forward.sam: Aligned forward reads
- alignment_reverse.sam: Aligned reverse reads



### 2. Unaligned reads extraction

The files generated in the previous step were used as inputs in order to extract unaligned reads with samtools,
then, bedtools bamtofastq is used for extracting FASTQ records from sequence alignments in BAM format:

```
samtools view -u -f 4 <alignment_forward.sam> > <alignment_forward.bam>
samtools view -u -f 4 <alignment_reverse.sam> > <alignment_reverse.bam>

bedtools bamtofastq -i <alignment_forward.bam> -fq <alignedReads_forward.fastq>
bedtools bamtofastq -i <alignment_reverse.bam> -fq <alignedReads_reverse.fastq>
```

### 3. Quality filtering with Prinseq


The unaligned reads in FASTQ format are filtered by quality with Prinseq, using the following parameters:

```
prinseq-lite.pl -fastq <alignedReads_forward.fastq> -trim_qual_right 25 -trim_qual_left 25 -trim_qual_type mean -trim_qual_window 5 -trim_qual_step 1  -trim_ns_left 1 -trim_ns_right 1 -ns_max_p 1 -ns_max_n 3 -min_qual_mean 25 -min_len 60 -out_good <clean_unal_forward.fastq> -out_bad <excluded_unal_forward.fastq>
prinseq-lite.pl -fastq <alignedReads_reverse.fastq> -trim_qual_right 25 -trim_qual_left 25 -trim_qual_type mean -trim_qual_window 5 -trim_qual_step 1 -trim_ns_left 1 -trim_ns_right 1 -ns_max_p 1 -ns_max_n 3 -min_qual_mean 25 -min_len 60 -out_good <clean_unal_reverse.fastq> -out_bad <excluded_unal_reverse.fastq>

```

This step will generate 4 files:

- clean_unal_forward.fastq and clean_unal_reverse.fastq: Forward and reverse clean reads
- excluded_unal_forward.fastq and excluded_unal_reverse.fastq: Forward and reverse excluded reads


### 4. Re-paired reads

Then, the pairs are joined and the orphan reads are discarded with the *in-house* script rePaired.py

```
python rePaired.py <clean_unal_forward.fastq> <clean_unal_reverse.fastq>
```

This step will generate 3 files:

- clean_unal_forward_paired.fastq: Forward paired reads
- clean_unal_reverse_paired.fastq: Reverse paired reads 
- clean_unal_forward_orphan.fastq: Orphan reads


At this point, it is advisable to run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) in order to guarantee the quality of the reads. 


### 5. Assembly with IDBA_UD


* First, the reads are joined and the pairs intercalated with the script called interleave_fq2fa.py:


The input files needed to run this script are:

- clean_unal_forward_paired.fastq: Forward paired reads
- clean_unal_reverse_paired.fastq: Reverse paired reads

    
```
python interleave_fq2fa.py <clean_unal_forward_paired.fastq> <clean_unal_reverse_paired.fastq> <clean_paired_interleaved.fa>
```


* Then, the assembly is performed with IDBA_UD using the interleaved file obtained in the previous step with the following parameters:

```
idba_ud -r <clean_paired_interleaved.fa> --mink 20 --maxk 120 --min_support 1 --min_contig 500 --pre_correction -o <clean_paired_interleaved_idba>


Allowed Options: 
  -o, --out arg (=out)                   output directory
  -r, --read arg                         fasta read file (<=128)
      --read_level_2 arg                 paired-end reads fasta for second level scaffolds
      --read_level_3 arg                 paired-end reads fasta for third level scaffolds
      --read_level_4 arg                 paired-end reads fasta for fourth level scaffolds
      --read_level_5 arg                 paired-end reads fasta for fifth level scaffolds
  -l, --long_read arg                    fasta long read file (>128)
      --mink arg (=20)                   minimum k value (<=124)
      --maxk arg (=100)                  maximum k value (<=124)
      --step arg (=20)                   increment of k-mer of each iteration
      --inner_mink arg (=10)             inner minimum k value
      --inner_step arg (=5)              inner increment of k-mer
      --prefix arg (=3)                  prefix length used to build sub k-mer table
      --min_count arg (=2)               minimum multiplicity for filtering k-mer when building the graph
      --min_support arg (=1)             minimum supoort in each iteration
      --num_threads arg (=0)             number of threads
      --seed_kmer arg (=30)              seed kmer size for alignment
      --min_contig arg (=200)            minimum size of contig
      --similar arg (=0.95)              similarity for alignment
      --max_mismatch arg (=3)            max mismatch of error correction
      --min_pairs arg (=3)               minimum number of pairs
      --no_bubble                        do not merge bubble
      --no_local                         do not use local assembly
      --no_coverage                      do not iterate on coverage
      --no_correct                       do not do correction
      --pre_correction                   perform pre-correction before assembly


```

* Finally, the script rename_contigs.py is used in order to change the contig names.

```
python rename_contigs.py <clean_paired_interleaved_idba.fa>
```


This step will generate a renamed FASTA file:

- clean_paired_interleaved_idba_renamed.fa


### 6. Contig analysis

The contigs are launched to the entire nucleotide NCBI database. Then, the best hit is obtained using blastDefaultBestHit_extract.py:


```
blastn -query clean_paired_interleaved_idba_renamed.fa -db nt -num_threads 12 -out clean_paired_interleaved_idba_renamed.fa.blastn
```

```
blastDefaultBestHit_extract.py -i clean_paired_interleaved_idba_renamed.fa.blastn


optional arguments:
  -h, --help  show this help message and exit
  -i BLAST    blast_file outfmt default

```


From this point forward, and with the information of all previous steps, the characterization of minicircles and maxicicles can be performed in order to
to try to circulate and find ORFs. In our study (see PAPER for more information), circularization was only performed for Minicircles and the ORFs were found in maxicicles contigs:


### 6.1 Minicircles circularization


The script minimusCircleFinder_v4.5.py is used in order to find circular contigs. This script use Minimus2 assembler internally. 


```
minimusCicleFinder_v4.5.py -i clean_paired_interleaved_idba_renamed.fa -k 120


optional arguments:
  -h, --help   show this help message and exit
  -i INPUT     Input file in fasta format (REQUIRED)
  -o OVERLAP   Minimum overlap, default 40
  -r REFCOUNT  Number of sequences is the first set, (Def 1)
  -c CONSERR   Maximum consensus error (0..1) (Def 0.06)
  -m MINID     Minimum overlap percentage identity for align. (Def 94)
  -t MAXTRIM   Maximum sequence trimming length (Def 20bp)
  -k KMER      max kmer used in the assembly (REQUIRED)

```

### 6.2 Maxicircle ORF-Finder

The possible ORFs of the maxicircle contigs are found using a script called leishORFinder_v1.py. The theoretical proteins are then analysed with blastp.


```
python leishORFinder_v1.py -i FASTA_FILE -l MIN_PROT_LEN


optional arguments:
  -h, --help       show this help message and exit
  -i FASTA_FILE    multifasta or single fasta file
  -l MIN_PROT_LEN  Minimun protein length

```


## MAINTAINERS

Current maintainers:
 * Alberto Rastrojo Lastras (CBMSO) - arastrojo@cbm.csic.es
 * Esther Camacho Cano (CBMSO) - ecamacho@cbm.csic.es 


## Authors

 * Esther Camacho Cano (CBMSO) - ecamacho@cbm.csic.es
 * Alberto Rastrojo Lastras (CBMSO) - arastrojo@cbm.csic.es
 * África Sanchiz Giraldo (CBMSO) - asanchiz@cbm.csic.es
 * Sandra González de la Fuente (CBMSO) - sandra.g@cbm.csic.es
 * Begoña Aguado Orea (CBMSO) - baguado@cbm.csic.es
 * José María Requena Rolanía (CBMSO) - jmrequena@cbm.csic.es

 
## CITE US

Please cite our paper published in [Genes](https://www.mdpi.com/journal/genes): Camacho, E., Rastrojo, A., Sanchiz, Á., Fuente, S. G.-D. L., Aguado, B., & Requena, J. M. (2019). [Mitochondrial Genomes: Maxicircle Structure and Heterogeneity of Minicircles](https://www.mdpi.com/2073-4425/10/10/758). Genes, 10(10), 758. doi: 10.3390/genes10100758 


## ACKNOWLEDGMENTS

This work was supported by grants (to B.A. and J.M.R.) from Proyecto del Ministerio de Economía, Industria y Competitividad SAF2017-86965-R, and by the Network of Tropical Diseases Research RICET (RD16/0027/0008); both grants are co-funded with FEDER funds. The CBMSO receives institutional grants from the Fundación Ramón Areces and from the Fundación Banco de Santander.




