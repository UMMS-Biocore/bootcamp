# Session 3 cont.: SAM/BAM alignment files

## Class Materials

You can follow the class materials below.

<b>1. Session 3.3: SAM/BAM Files</b><br />

<div align="left">
  <a href="https://www.youtube.com/watch?v=56BIJEntE8w"><img src="https://img.youtube.com/vi/56BIJEntE8w/0.jpg" alt="Session 3.3"></a>
</div>

Slides: [SAM/BAM file formats](Session3.3_general_pipeline_alignment_files.pptx)

## SAM files

SAM is an acronym for Sequence Alignment/Mapping

It is not common to find a sam file since they take up a lot of space, but here is a small example. Since it is a text file, you can view it just as you had been viewing all other text files we have worked with:

```
$ head /project/umw_biocore/class/mm_example.sam
```

## BAM files

BAM files are binary, and for that reason, they can't be easily manipulated as other text files. The goal of this tutorial is to introduce you to the tools needed to manipulate them.

### SAMTools commands

In order to work with BAM files, you can use the SAMTools library of utilities. Below is a link to the full description of the commands:

http://www.htslib.org/doc/samtools.html

This is already installed in our cluster. To use it, **load the module for the most recent version** (remember you can query all the modules available to find it), and follow the exercises below:

- Get the help for samtools commands:

```
$ samtools
```

Notice that there are multiple commands to use and they require you to first type samtools, followed by the subcommand:

```
#example
samtools view
samtools sort
samtools merge
```

- Use the view command to print the first 10 alignments of this bam file:

```
$ samtools view /project/umw_biocore/class/mm_example.bam | head
```

- Count the number of records in this file:

```
$ samtools view /project/umw_biocore/class/mm_example.bam | wc -l
```

- View the header only for the same bam file (by default only reads are returned)

```
$ samtools view -H /project/umw_biocore/class/mm_example.bam | head
```

Headers may or may not be present in SAM/BAM files, and they always start with the @ character. The last line in the header includes the command used to generate the bam file and is useful to store a record of the parameters used for the alignment!

### Homework (Advanced)

- Using a single command line (may include multiple commands using **pipes**), print the last read in the bam file used above in a FASTQ format. Save a screenshot of the command and result.

_hint: awk allows you to change the output separator of fields or print any string you want, and a newline can be specified by_ "\n"

```
$ head -1 ~/bootcamp/RNA-Seq/reads/exper_rep1.1.fq | awk '{print "Hello World\n"$1}'
```

Expected Result:

```
@NS500602:410:HJGVMBGX2:2:23305:25539:14380
TTTGCATCGGTGTAATAGGGCACCAGTGACTCAGGGGGAAGC
+
EEEEEEEEEEEEEEEEEEEEEEAAEEEEEEEAEEEEEAAAA/
```
