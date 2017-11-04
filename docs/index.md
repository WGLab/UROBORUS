## Introduction

UROBORUS is a computational pipeline to detect circular RNA from RNA-Seq data, based on junction reads from back spliced exons.

## Usage

Before using UROBORUS.pl, you should use TopHat to align the reads to genome, and get the unmapped.sam file.

```
usage:
perl UROBORUS.pl -index /path/genome -gtf /path/genes.gtf -fasta /path unmapped.sam accepted_hits.bam
Options:
-index: genome index (use bowtie1 index);
-gtf:   gene annotation file (*.gtf file);
-fasta: path for genome sequence in fasta file (*.fa) in separate chromosome;
-p:     threads (Integer, default = 6);
-temp:  keeping the temporary file;
-help:  usage help;
```

## Note

1. If the genome sequence in each chromosome (chr1.fa, chr2.fa, chr3.fa)is saved in the directory `/home/circRNA/Gene`, the path for genome sequence should be set as `-fasta /home/circRNA/Gene`;
2. The following files ( genome index, gene annotation, genome sequence) should be download from TopHat webstie: http://ccb.jhu.edu/software/tophat/igenomes.shtml


## Example

```
perl UROBORUS.pl -index /home/***/index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome -gtf /home/***/circRNA/Gene/genes.gtf -fasta /home/***/circRNA/Gene unmapped.sam accepted_hits.bam
```

## Prerequisites

Software Prerequisites:

The following three software should be installed in your cluster or computer before running the UROBORUS.pl

* TopHat

`tophat -p 6 -o RL_9_tophat_out /home/***/index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome
/home/***/data/Glioblastoma/RL_9_OLI_GTCCGC_L008_R1_001.fastq /home/***/data/Glioblastoma/RL_9_OLI_GTCCGC_L008_R2_001.fastq`

* samtools
Convert unmapped.bam to unmapped.sam (using samtools view)

`samtools view unmapped.bam > unmapped.sam`

* Bowtie1
Within UROBORUS.pl, Bowtie1 should be used.
However, before running UROBORUS.pl, You run Tophat to get unmapped.sam file, Bowtie1 or Bowtie2 should be selected based on read length. 

## Output file format

The first 3 columns are the same with bed file format.

1. Chromosome
2. start of junction
3. end of junction
4. strand
5. Parental gene name
6. genomic distance
7. read counts
8. matched transcript id

## Reference

Song X, Zhang N, Han P, Lai RK, Wang K, Lu W. Circular RNA Profile in Gliomas Revealed by Identification Tool UROBORUS. Nucleic Acids Research, 2016, 44:e87.

## Contact

Please contact Xiaofeng Song (xfsong@nuaa.edu.cn) for questions and comments.

<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-73567476-1', 'auto');
  ga('send', 'pageview');

</script>
