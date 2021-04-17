Session 2 cont.: Usefull commands and tools
========

Expected learning outcome
========

We will learn some usefull commands that will help you to use linux. We will also use some editors to create some text files. These editors can help you to manupilate text files and create bash scripts to execute commands in more organized way.

Overview
========
  * [Class Materials](#class-materials)
  * [Introduction](#introduction)
  * [Processes Management](#processes-management)
  * [Pipes](#pipes)
  * [grep Command](#grep-command)
  * [AWK Command](#awk-command)
  * [vi Editor](#vi-editor)
  * [Session Homework](#session-homework)

## Class Materials
You can follow the class materials below.

<b>1. Session 2.2: Useful linux commands and tools</b><br />

<div align="left">
  <a href="https://www.youtube.com/watch?v=vrwklLU-Ic8"><img src="https://img.youtube.com/vi/vrwklLU-Ic8/0.jpg" alt="Session 2.2"></a>
</div>

## Introduction
There are a couple of more commands that is required to learn to be more comfortable in linux environment. In this session we will learn processes in management, pipe usage, grep, awk, sed commands that are usefull for text manupilation and text search in files and vi editor.

## Processes Management
When you run a command, it is called a process and there are two ways you can run it −

- Foreground Processes
- Background Processes

### Foreground Processes
By default, every process that you start, runs in the foreground. For example, when you start "sleep 300" command. It will wait 300 seconds and you cannot do anything else while it is running. 
	
	$ sleep 300
		
### Background Processes
To send a process to the background, you can start the command and use Ctrl+Z keys to send it to background. When you press these keys, it will stop the run. You need to run `bg` command to continue to the run.

	$ sleep 300
	
press Ctrl+Z in your keyboard

	$ bg

This command will continue running the command.

Or you can run the process in the background with `&` argument.

	$ sleep 300 &

This will send the process to the background without Ctrl+Z and bg command.

If you want to get the process, running in the background, to the foreground use `fg` command.

	$ fg 

* **In order to exit and stop the command use "Ctrl+C"**

### Listing Running Processes
To see running processes you can use `ps` command or `ps -f` to show it with full details.

	$ ps

or
	
	$ps -f
	UID   PID  PPID   C STIME   TTY           TIME CMD
	501  1384  1383   0  2:29AM ttys004    0:00.23 -bash
	501  4373  1384   0  1:56PM ttys004    0:00.00 sleep 300
	501  1395  1394   0  2:29AM ttys005    0:00.02 -bash
  
| Column | Description                                                  |
| ------ | ------------------------------------------------------------ |
| UID    | User ID that this process belongs to (the person running it) |
| PID    | Process ID                                                   |
| PPID   | Parent process ID (the ID of the process that started it)    |
| C      | CPU utilization of process                                   |
| STIME  | Process start time                                           |
| TTY    | Terminal type associated with the process                    |
| TIME   | CPU time taken by the process                                |
| CMD    | The command that started this process                        |

### Stopping Processes
If the process is running in the foreground you can use Ctrl+C command to stop it or you can use `ps` command to learn PID of that process and use `kill` command to stop it. For example to kill `sleep 300` command above, first check the PID which is reported at above as 4373. The command below will kill it.

	$ kill 4373
	Terminated
	
If a process ignores a regular kill command, you can use `kill -9` followed by the process ID as follows.

	$ kill -9 4373
	Terminated
	
* **Note that, PID or a Process ID is different than JOBID in LSF schedular.**

### The top Command
The top command is a very useful tool for quickly showing processes sorted by various criteria.

It is an interactive diagnostic tool that updates frequently and shows information about physical and virtual memory, CPU usage, load averages, and your busy processes.

	$ top

* You can exit from the "top" command by typing "Ctrl+C"

## Pipes

You can connect two commands together so that the output from one program becomes the input of the next program. Two or more commands connected in this way form a pipe.

To make a pipe, put a vertical bar `|` on the command line between two commands.

If you want to count the # of lines of an output of `ls -l` command, you can use the command below

	$ ls -l | wc -l

Or if you want to count # of lines in a fastq file that we copied into the directory below;
These were paired end files, so the # of lines in between control_rep1.1.fq and control_rep1.2.fq has to be the same. You can check that using the commands below;

	$ cd ~/bootcamp/RNA-Seq/reads/
	$ cat control_rep1.1.fq | wc -l
	99152
	$ cat control_rep1.2.fq | wc -l
	99152

If they aren't equal, you might need to re-download the files. 

## The grep Command

The grep command searches a file or files for lines that have a certain pattern. 

	$ grep pattern file(s)

The simplest use of grep is to look for a pattern consisting of a single word.

	$ ls -l > ~/mylist.txt
	$ grep Apr ~/mylist.txt

The same command could be run with pipe without writing the outout of the "ls -l" command to  a file.

	$ ls -l | grep Apr 

| Option | Description                                                     |
| ------ | --------------------------------------------------------------- |
| -v     | Prints all lines that do not match pattern.                     |
| -n     | Prints the matched line and its line number.                    |
| -l     | Prints only the names of files with matching lines (letter "l") |
| -c     | Prints only the count of matching lines.                        |
| -i     | Case insensitive search. Matches either upper or lowercase.     |

Let's play with grep on our sequence files.

<img src="images/grep.png">

If it takes more than one page you can use less command.

	$ cd ~/bootcamp/RNA-Seq/reads
	$ grep AAA --color=always control_rep1.1.fq | more

* Tip: Press `q` to exit from "more".

You can also use `--color` in short. Let's search another sequence CAGAGTTC and put the results into a separate file.

	$ cd ~/bootcamp/RNA-Seq/reads
	$ grep CAGAGTTC --color control_rep1.1.fq > ~/bootcamp/seq.txt

To learn how many reads include CAGAGTTC in control_rep1.1.fq file

You could use "-c" argument;
	
	$ grep CAGAGTTC -c control_rep1.1.fq
	
or with pipe and `wc` command:

	$ grep CAGAGTTC control_rep1.1.fq | wc -l

Not every command has -c argument, so learning second method might be useful for other commands like `cat`.

If you want to search a motif in your reads starts with CAGAGT and ends with TGA. There could be anything in between. To represent any nucleotide in between these two sequences, you can use `.*`. So the pattern will be `CAGAGT.*TGA`. 
 
	$ grep CAGAGT.*TGA --color control_rep1.1.fq

If you want to eliminate the reads that include the pattern above but include another sequence
GCT...CGCG. "..." means any three nucleotides.	 
	
	$ grep -v CAGAGT.*TGA control_rep1.1.fq | grep GCT...CGCG --color

To learn more about it

<https://ryanstutorials.net/linuxtutorial/grep.php>

Cheat sheet:

<https://ryanstutorials.net/linuxtutorial/cheatsheetvi.php>

## AWK command

**awk** is an interpreted language designed to process text files. **awk** is especially good for the files and command outputs that have columns. You can access any column of a file easily with awk.

<img src="images/awk.png">

| Content     | Awk variable |
| ----------- | ------------ |
| Entire Line | $0           |
| Column 1    | $1           |
| Column 2    | $2           |
| .           |              |
| .           |              |
| Column i    | $i           |
| Line Number | NR           |

For example, if we only want to list username and the name of the directories we can use these commands;

	$ cd ~/bootcamp/RNA-Seq
	$ ls -l
	total 38
	drwxrwxr-x 2 ak97w umw_manuel_garber  51 Apr  7 23:44 mm10
	drwxrwxr-x 2 ak97w umw_manuel_garber  34 Apr  7 23:46 progs
	drwxrwxr-x 2 ak97w umw_manuel_garber 408 Apr  7 23:32 reads

Let's list column 3 which is user id and column 9 and put space in between these columns.

	$ ls -l > mylist.txt
	$ awk '{print $3" "$9}' mylist.txt
	
Or without writing the list into a file, we could redirect the output of "ls -l" directly to the awk command.
	
	$ ls -l | awk '{print $3" "$9}'

To put a space between columns, I used space in quotes `" "`. If you want to put a **tab** character in between columns you can use `"\t"` this usage is called using "escape character" that you might see it in other tutorials;

	$ ls -l | awk '{print $3"\t"$9}' 

### Getting all nucleotide sequences from a fastq file

In fastq files, there are 4 lines for each read. The nucleotide sequence of the reads is on the second line for every forth lines. We can get them using a very simple modular arithmetic operation,

	$ cd ~/bootcamp/RNA-Seq/reads/
	$ head -n 8 control_rep1.1.fq 
	
	@HWI-ST333_0273_FC:7:1205:9843:15143#CCGAGT/1
	CAAGGAAGCACATGACCGAGCAGAAAATACCCAGTTTGTC
	+
	iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
	@HWI-ST333_0273_FC:7:2116:8217:37625#CCGAGT/1
	GTCAGCTTCCTGATGTTCTCCAGGCCACTGTACACTACAT
	+
	iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

In order to get the nucleotide sequence of the reads we can use awk command. `NR` here is the line number `NR % 4 == 2` will use the modulus operator to get second lines of every forth lines in a file.  

	$ awk '{if(NR % 4 == 2)print($0)}' control_rep1.1.fq | less
	
### An excercise to kill multiple jobs in the cluster using awk

Let's submit three dummy jobs

	$ bsub -J job1 "sleep 300"
	$ bsub -J job2 "sleep 300"
	$ bsub -J job3 "sleep 300"

If you need more time than 300 seconds to finish this section, you can increase it.
Now we can build the necessary command to kill those dummy jobs.

	$ bjobs 
	JOBID      USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
	4358754    ak97w   RUN   short      ghpcc06     c41b01      job1       Apr 14 17:36
	4358755    ak97w   RUN   short      ghpcc06     c41b01      job2       Apr 14 17:36
	4358756    ak97w   RUN   short      ghpcc06     c41b01      job3       Apr 14 17:36

To print the jobs where their job names include "job", we can use a command below;

	$ bjobs | grep job
	4358754    ak97w   RUN   short      ghpcc06     c41b01      job1       Apr 14 17:36
	4358755    ak97w   RUN   short      ghpcc06     c41b01      job2       Apr 14 17:36
	4358756    ak97w   RUN   short      ghpcc06     c41b01      job3       Apr 14 17:36

As you can see, it skipped the first line (header line). Now we can print the job IDs.

	$ bjobs | grep job| awk '{print $1}'
	4358754
	4358755
	4358756

To prepare kill commands for each job;

	$ bjobs | grep job| awk '{print "bkill "$1}'
	bkill 4358754
	bkill 4358755
	bkill 4358756

The commands above will only write the bkill commands to the screen but it won't execute. To execute we need to send the output to bash like below;

   	$ bjobs | grep job| awk '{print "bkill "$1}'|bash
	Job <4358754> is being terminated
	Job <4358755> is being terminated
	Job <4358756> is being terminated
	
This is one of the many ways to do the same operation. 

You can check other tutorials and cheat sheets online;

<https://www.grymoire.com/Unix/Awk.html>

<http://www.awklang.org/asset/krumnisCheatSheet.pdf>

## The vi Editor 

There are many ways to edit text files in linux. **vi** is one of the most powerful editor you can use from command line. To start the **vi** editor you use 
	
	$ vi testfile
	
The above command will generate the following output;

	|
	~
	~
	~
	~
	~
	~
	~
	~
	~
	~
	~
	~
	"testfile" [New File]    
 
 A tilde represents an unused line.

### Operation Modes

- ***Command mode:*** As you enter vi, command mode will be active. This mode enables you to perform administrative tasks such as saving the files, executing the commands, moving the cursor, cutting (yanking) and pasting the lines or words, as well as finding and replacing. In this mode, whatever you type is interpreted as a command.
- ***Insert mode:*** If you press "i" in the **command mode**, it will switch to **insert mode**. This mode enables you to insert text into the file. Everything that's typed in this mode is interpreted as input and placed in the file. Press the **Esc** key twice to return back to command mode after you complete the editing.

***Hint:*** If you are not sure which mode you are in, press the Esc key twice; this will take you to the command mode.

### Getting Out of vi

First you need to return to **command mode** by press the **Esc** key twice. 
* If you haven't edited the file, type `:q` (press enter to exit)
* If you've edited the text and want to quit without saving the changes, type `:q!` (press enter to exit)
* Save and exit is `:wq` (press enter to exit)
* Only save the text is `:w`

### Editing Files
To edit the file, you need to be in the insert mode. There are many ways to get in to insert mode but the easiest way to get into insert mode please press i. To exit from insert mode press **ESC**. 

When you are in instert mode you will see -- INSERT -- at the bottom of the page;

	|
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	~                                                                                                                                                      
	-- INSERT --

You can start writing now and to save and exit. Press ESC and :wq

Let's write a command `ls -l` into a file with vi.

	$ cd ~
	$ vi myls
 
 Go into insert mode write `ls -l` into this file and save and exit (Press ESC and :wq).
 
 Lets make this file executable.
 
 	$ chmod +x myls
 	
 Now it is ready to execute.
 
 	$ ./myls
 	
Or use the **full path** where the file is located. If this file is created in your home

	$ ~/myls
	
It will run "ls -l" command.

* Tip: If your myls file located somewhere else, you can learn the **full path** of a directory by running `pwd` command.
 
There are many commands to search, replace, delete the lines etc. in vi. To learn more about vi, you can use the tutorial below;

<https://ryanstutorials.net/linuxtutorial/vi.php>
 
* **While using vi, make sure the file is a small text file, bash script or a program. Don't open big files (e.g. fastq) with vi. It can open but takes a lot of time and not advised**


## Session Homework:

1. Submit two dummy jobs to the long queue and three dummy jobs to the short queue. Then get a list of your jobs that have been submitted to the long queue only. 

	- Hint 1: Use sleep 300 to create a dummy job.
	- Hint 2: Use bsub to submit a job. Remember that -q parameter is used to specify the queue.
	- Hint 3: Recall that bjobs can be used to list your jobs in the cluster.
	- Hint 4: Use what you have learned so far and put the pieces together.

2. Print the "start codon" location of Fgf21 gene using the annotation file(gtf) in the following format and send the command to one of us.
	
	The annotation file location:
	
		/share/data/umw_biocore/genome_data/mouse/mm10/ucsc.gtf 
	
	The output format:
		
		chr7	"Fgf21";	start_codon	45615304
		
	* Hint-1: Use grep and awk commands with pipe.
	* Hint-2: You can find "start codon" by searching `start_codon`.

 
