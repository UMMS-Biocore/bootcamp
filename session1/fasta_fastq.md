# Session 3 cont.: FASTA and FASTQ files

## Class Materials

You can follow the class materials below.

<b>1. Session 3.2: Genomic Resources</b><br />

<div align="left">
  <a href="https://www.youtube.com/watch?v=gc-wIrxpSMw"><img src="https://img.youtube.com/vi/gc-wIrxpSMw/0.jpg" alt="Session 3.2"></a>
</div>

Slides: [Fasta/Fastq file formats](Session3.2_sequence_file_formats.pptx)

## FASTA files

Now you learned the location of the genomes in the cluster, find a fasta file for the mouse reference (mm10)

```
$ ls -l /share/data/umw_biocore/genome_data/mouse/mm10/*.fa
-rwxr-xr-x 1 gxy11w umw_manuel_garber 2785490220 Apr 29  2014 /share/data/umw_biocore/genome_data/mouse/mm10/mm10.fa
-rwxr-xr-x 1 gxy11w umw_manuel_garber   99845877 Aug  2  2019 /share/data/umw_biocore/genome_data/mouse/mm10/rsem_ref.idx.fa
-rw-r--r-- 1 gxy11w umw_manuel_garber   95664752 Jul 10  2019 /share/data/umw_biocore/genome_data/mouse/mm10/rsem_ref.n2g.idx.fa
-rwxr-xr-x 1 gxy11w umw_manuel_garber   99845877 Aug  2  2019 /share/data/umw_biocore/genome_data/mouse/mm10/rsem_ref.transcripts.fa
```

Now look at the first 10 lines of the mouse mm10.fa file using **head**

```
$ head /share/data/umw_biocore/genome_data/mouse/mm10/mm10.fa
```

And the last 10 lines of the same file using **tail**

```
$ tail /share/data/umw_biocore/genome_data/mouse/mm10/mm10.fa
```

Looking at the output, notice two more details about FASTA files:

- They can contain the N character in nucleotide sequences, not just A/C/T/G. Usually this is used to represent gaps and unknown sequences in the genome assembly, for example **telomeres** and **centromeres**.
- Sometimes the sequence will be encoded as lowercase actg in place of the standard ACTG notation, for most uses this will not make a difference, but some programs can interpret these as different information, and it is commonly used to note the regions of the genome with **repetitive elements** so that they can be easily **_masked_** if that is desired.

## FASTQ files

During the homework exercise in [Session 1: Linux environment](session1/session1.md), you made a copy of several fastq files to your home directory:

```
$ ls -l ~/bootcamp/RNA-Seq/reads/
total 32946
-rwxr-x--- 1 ed70w umw_manuel_garber 3233830 Apr 15 15:45 control_rep1.1.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 3233830 Apr 15 15:45 control_rep1.2.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2112249 Apr 15 15:45 control_rep2.1.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2112249 Apr 15 15:45 control_rep2.2.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2283880 Apr 15 15:45 control_rep3.1.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2283880 Apr 15 15:45 control_rep3.2.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2356860 Apr 15 15:45 exper_rep1.1.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2356860 Apr 15 15:45 exper_rep1.2.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2807338 Apr 15 15:45 exper_rep2.1.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 2807338 Apr 15 15:45 exper_rep2.2.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 1251985 Apr 15 15:45 exper_rep3.1.fq
-rwxr-x--- 1 ed70w umw_manuel_garber 1251985 Apr 15 15:45 exper_rep3.2.fq
```

#### Exercise

- Take a look at the first and last 8 lines in one of these files

- Count how many reads are in the file **exper_rep1.2.fq**. Since there are 4 lines per read in this fastq, you can divide the result of the number of lines by 4 using awk to get the correct number:
  ```
  $ wc -l exper_rep1.2.fq | awk '{print $0/4}'
  ```

### Homework

During [Session 2 Extra: Useful commands and tools](../session2/usefull.md) you became familiar with the command **grep**, but it has extra power when used with certain arguments:

```
$ grep --help
```

- What is the nucleotide sequence of the read named HWI-ST333_0273_FC:8:1212:9473:35651#CTGGGC/1 in the file **exper_rep2.1.fq**

Expected result:

```
AAAGCGGTCAACTGGAAACTAAATGGAGCAGCAGCTCCTC
```
