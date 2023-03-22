# Session 3: Genomic resources in the cluster

## Class Materials

You can follow the class materials below.

<b>1. Session 1.3: Genomic Resources</b><br />

Slides: [Introduction to Sequencing Analysis](Session1_sequencing_introduction.pptx)

## Finding genomes

All the reference files for **genomes** and **annotations** hosted by the UMMS Biocore can be found here:

```
/share/data/umw_biocore/genome_data/
```

List all the directories, notice how there are entries for each organism

```
$ ls -l /share/data/umw_biocore/genome_data/
```

In each organism directory, there can be one or more reference genome assemblies

```
$ ls -l /share/data/umw_biocore/genome_data/mouse/
total 192
drwxrwxr-x 4 onur.yukselen-umw alper.kucukural-umw 12288 Sep 18  2020 mm10
drwxrwxr-x 3 onur.yukselen-umw alper.kucukural-umw  8192 Jul 24  2020 mm10_gencode_m25
drwxrwxr-x 3 onur.yukselen-umw alper.kucukural-umw  4096 Nov  5  2019 mm10_meta
drwxrwxr-x 3 onur.yukselen-umw alper.kucukural-umw  8192 Apr 25  2019 mm9
```

Using the human (hg38) reference, we can easily see how many chromosomes are in this annotation and their length in nucleotides:

```
$ cat /share/data/umw_biocore/genome_data/human/hg38_gencode_v28/chrNameLength.txt
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
chr11	135086622
chr12	133275309
chr13	114364328
chr14	107043718
chr15	101991189
chr16	90338345
chr17	83257441
chr18	80373285
chr19	58617616
chr20	64444167
chr21	46709983
chr22	50818468
chrX	156040895
chrY	57227415
chrM	16569
```

#### Exercise:

Count the chromosomes (file with the same name) using the most recent annotations for _D. melanogaster_ and rat

_hint: one line per chromosome in the file_
