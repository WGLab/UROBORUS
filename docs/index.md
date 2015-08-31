# UROBORUS

UROBORUS is a computational pipeline to detect circular RNA supported by junction reads from back spliced exons.

## Download

Download the latest version of UROBORUS. Version: 1.0.0. Last Modified: 2015-8-05
Download here

## Usage

Authors: Xiaofeng Song (xfsong.nuaa@gmail.com or xfsong@nuaa.edu.cn)
Usage:
Before using UROBORUS.pl, you should use TopHat to align the reads to genome, and get the unmapped.sam file.
UROBORUS.pl 1.0.0 #circRNA identification tool in total RNA-seq data
usage:
perl UROBORUS.pl -index /path/genome -gtf /path/genes.gtf -fasta /path unmapped.sam
Options:
-index:	genome index (use bowtie index);
-gtf:	gene annotation file (*.gtf file);
-fasta:	path for genome sequence in fasta file (*.fa) in separate chromosome;
-p:	threads (Integer, default = 6);
-temp:	keeping the temporary file;
-help:	usage help;

## Note:

1.If the genome sequence in each chromosome (chr1.fa, chr2.fa, chr3.fa, ¡­)is saved in the directory /home/circRNA/Gene, the path for genome sequence should be set as ¡°-fasta /home/circRNA/Gene¡±;
2.The following files ( genome index, gene annotation, genome sequence) should be download from TopHat webstie: http://ccb.jhu.edu/software/tophat/igenomes.shtml


## Example:

```
perl UROBORUS.pl -index /home/***/index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome -gtf /home/***/circRNA/Gene/genes.gtf -fasta /home/***/circRNA/Gene unmapped.sam
```

## Installation

Software Prerequisites:

The following three software should be installed in your cluster or computer before running the UROBORUS.pl

* TopHat

```
tophat -p 6 -o RL_9_tophat_out /home/***/index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome
/home/***/data/Glioblastoma/RL_9_OLI_GTCCGC_L008_R1_001.fastq /home/***/data/Glioblastoma/RL_9_OLI_GTCCGC_L008_R2_001.fastq
```

* samtools

Convert unmapped.bam to unmapped.sam (using samtools view)

```
samtools view unmapped.bam > unmapped.sam
```

* Bowtie

## Output file format:

The first 3 columns are the same with bed file format.
1. Chromosome
2. start of junction
3. end of junction
4. strand
5. Parental gene name
6. genomic distance
7. read counts
8. matched transcript id


