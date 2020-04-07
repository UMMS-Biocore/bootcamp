Session 1: Linux/Unix Environment
========

Disclaimer
-------

For the sake of clarity, the concepts mentioned in this tutorial have been simplified significantly. Linux is not Unix, but it is a Unix-like operating system. Linux system is derived from Unix. Similar commands you will learn in this tutorial may be applicable to other unix systems.
Most of the command line tools and programs we explain have many more features that cannot be found here. Linux command line tools are very well documented and their manual pages are freely available.
These notes are prepared for the users of UMMS Cluster. UMMS Cluster system has features that might not exist in standard Unix/Linux systems such as loading modules.
A rigorous treatment of topics on Unix and Bash can be found in various books in several levels.

Expected learning outcome
========

To understand the basics of Linux environment, how to use command line and to familiarize yourself with some standard tools and commands using a terminal.

Overview
========
  * [Introduction](#introduction)
  * [Before you start](#before-you-start)
  * [Getting Started](#getting-started)
  * [Basic Commands](#basic-commands)
 

## Introduction

Linux is an operation system like Microsoft Windows or OSX. An operating system is a collection of software that manages system resources (such as memory, cpu(s), display and etc.) and provide applications a simpler interface to the system
hardware.

<img src="images/operating_system.png">
<img src="images/unix.png">
<img src="images/unix_like.png">
<img src="images/unix_family.png">
<img src="images/terminal.png">

### Unix shell

A shell is a software that runs inside a terminal and interprets and executes user commands.
One of the most popular shells being used today is called BASH (Bourne Again Shell). In this tutorial we will use bash. There are other shells like zsh or sh. They can also be used if you need for a specific functionality that is supported by that shell. 

## Before you start
### Convention: 

***$*** denotes a prompt for the command line. It is not meant to be typed while running the commands in this tutorial. All Bash commands will appear in a gray box that you will run in your terminal. In this example below, you will only write "your\_command" in the terminal an press enter. It will then execute the command in your terminal. (Please only run commands in the lines starts with \"$\" symbol.)

	$ your_command

Let's start running some commands. First command will be "ssh" to connect to the UMASS cluster.

	$ ssh username@ghpcc06.umassrc.org
	

Let’s verify that we are at the right place. In this case only run "hostname", not the output of the command (ghpcc06).

	$ hostname
	ghpcc06
	
You need to see "ghpcc06" in your terminal as an output. ghpcc06 is our "head node". You will learn what head or child nodes are in these tutorials later. 


## Getting started

To print something on the screen. We use echo command. 
	
	$ echo Hello World
	
If you just run 

	$ Hellow World
	-bash: Hello: command not found
	
You will get this error "command not found". Because there is no any command called "Hello". As you can see, this command is interpreted by "bash" and when there is an error, -bash reported the error.


## Basic Commands
  * [List a directory](#list-a-directory)
  * [Locating Applications and Software Packages](locating-applications-and-software-packages)
  * [Modular system](#modular-system)
  * [Print working directory](print-working-directory)
  * [Accessing Directories](accessing-directories)
  * [Exploring the Filesystem](exploring-the-filesystem)
  * [Creating a directory](creating-a-directory)

Before we explore the commands used to manipulate the Linux environment, we should take a quick look at the structure of the environment itself.

You can think of a Linux file system as an upside-down tree. See the diagram below. In this diagram, the boxes are directories or folders. "script.sh" is a file.
 
<img src="images/dir_structure.png">

At top you will see the following symbol; "/". It is called "root" directory. For example if user1 wants to access Documents directory,  user will use "/home/user1/Documents" when it is needed. This directory structure is called "Full Path" of a directory. Full path of script.sh in this example will be "/home/user1/script.sh". 

### List a directory
Probably the most often used command in Linux is the "ls" command. It is used to list the contents of a directory. 
Unlike many other operating systems, Linux is case-sensitive. In other words, if you type "LS" instead of "ls", Linux will not recognize the command. This applies to director and file names, like "home" and "script.sh", as well.

Lets run some commands below one at a time and list the content of the directory. 

	$ ls /
	$ ls /project
	$ ls /nl
	$ ls /project/umw_biocore/class/

For example, last command result will be

    $ ls /project/umw_biocore/class/
    class.html  class.R  data.tsv  funcs.R

We will visit ls command later a couple of more time to learn about the size of a file or permissions, creation time, or owner of file or a directory.

### Locating Applications and Software Packages
Depending on the specific distribution, programs and software packages can be installed in various directories. In general, executable programs should live in the following directories

```
/bin
/usr/bin
/sbin
/usr/sbin
/opt.
```

One way to locate programs is to employ the ``which`` utility. For example, to find out exactly where the diff program resides on the filesystem:

	$ which diff
	/usr/bin/diff

If which does not find the program, whereis is a good alternative because it looks for packages in a broader range of system directories:

	$ whereis diff
	diff: /usr/bin/diff /usr/share/man/man1/diff.1.gz

### Modular system

Other software packages we will use in our cluster is installed in a modular structure to support multiple versions of a software in a single system. For example, "STAR" is a program  for splice aware genomic alignments. However, none of the STAR versions is loaded in the system. 

	$ which STAR
	
	/usr/bin/which: no STAR in (...)

which command couldn't find the STAR command in any of our paths.

Please write the command below and press "tab button" in your keyboard. It will list all available version of star. (Tab can also be used for auto-completion of any other command.)

	$ module load star/2.   #(Press tab button here)
	star/2.3.0e  star/2.4.2a  star/2.5.3a  star/2.7.0e  

To load the desired module  

	$ module load star/2.7.0e  

Run "which" command again to see the location of the command;

	$ which STAR
	
To unload a module;

	$ module unload star/2.7.0e  

To list all available modules in the cluster;

	$ module avail 

If a version of a software package you want to use is not installed in our cluster. You can always send an email to hpcc-admins to have it installed. 


### Print working directory

pwd command will show what directory you are in.

To see the working directory,

	$ pwd


### Accessing Directories
The following commands are useful for directory navigation:

|Command|Result|
|-------|-----------|
|cd 	|Change to your home directory|
|cd ~	|Change to your home directory (same like above)|
|cd ..|Change to parent directory|
|cd - |Change to previous directory|
|cd /	|Changes your current directory to the root (/) directory|
|cd /project/umw_biocore/class|Changes your current directory to a specific directory|

Please use all the commands in the table above and run "pwd" command after that to see which directory you are in and "ls" to see the content of the directory.

	$ cd
	$ pwd
	$ cd ..
	$ pwd
	$ cd -
	$ pwd
	$ ls
	$ cd /
	$ pwd

### Exploring the Filesystem
The tree command is a good way to get a bird’s-eye view of the filesystem tree. The following commands can help in exploring the filesystem:

|Command|Result|
|-------|-----------|
|ls 	  |List the contents of the present working directory|
|ls –a  |List all files including hidden files and directories|
|ls –l  |Detailed list of files and directories|
|ls –lh  |Detailed list of files and directories where the file sizes are reported in human readable format|
|tree   |Displays a tree view of the filesystem|
|tree -d|Just list the directories and suppress listing file names|


Please use the commands in the table using the example directory "/project/umw_biocore/class". You can use these command with any other directories too. 

To see the complete list of options of a command. 

	$ man ls
	
***To exit from manual of a command, please press "q" button in your keyboard.***

Try commands below;

	$ ls /project/umw_biocore/class
	$ ls –a /project/umw_biocore/class

Which file is different in the output of both commands above?

The hidden files start with ".". 

	$ ls –l /project/umw_biocore/class
	$ ls –lh /project/umw_biocore/class
    
You can also run -l -a command together
    
	$ ls -a -l /project/umw_biocore/class
    
ls command can also recognize the parameters together like in -lh suports in any order
 
 	$ ls -al /project/umw_biocore/class
 	
tree command shows all files and folders recursively. If you have many files just press Ctrl+C, in the directory you want to check the tree.

	$ tree /project/umw_biocore/class

***If any of your commands run too long. You can always use Ctrl+C in your keyboard to stop the execution.***
	
List the Tree with only the folder names

	$ tree -d /project/umw_biocore/class
	
### Creating a directory

Let's first go to our home directory.

	$ cd
	
Create a directory called "bootcamp"

	$ mkdir bootcamp
	
Lets go into this directory and print where we are.

	$ cd bootcamp
	$ pwd

You can try going to another directory using the full path

	$ cd /project/umw_biocore/class
	$ pwd
	$ ls
	
Please go back using full path, or using "~" symbol. "~" means your home directory.


	$ cd ~/bootcamp
	$ pwd
	
Or you could have used 

	$ cd /home/your_user_id/bootcamp
	
The directory could also be created like below. 
	
	$ mkdir ~/bootcamp

Or using the full path

	$ mkdir /home/your_user_id/bootcamp

However, if you try the commands above, bash will give an error, since directory is already created.


If you want to create a directory inside of another directory which is not created before, you can use "-p" command. For example, creating ~/bootcamp/dir1/dir2/dir3

	$ mkdir ~/bootcamp/dir1/dir2/dir3
	mkdir: cannot create directory `/home/your_user/bootcamp/dir1/dir2/dir3': No such file or directory

This command will give an error.

When you use -p option will create all those directories;

	$ mkdir -p ~/bootcamp/dir1/dir2/dir3
	$ cd ~/bootcamp/dir1/dir2/dir3
	$ pwd
	
You can also create multiple directories with one command.

	$ mkdir ~/bootcamp/first ~/bootcamp/second
	
If you go to that directory and get a list;

	$ cd ~/bootcamp
	$ ls
	dir1  first  second
	
You can also use tree command

	$ tree ~/bootcamp
	/home/your_user/bootcamp
	├── dir1
	│   └── dir2
	│       └── dir3
	├── first
	└── second

