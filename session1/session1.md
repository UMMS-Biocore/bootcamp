Session 1: Linux Environment
========

Disclaimer
-------

For the sake of clarity, the concepts mentioned in this tutorial have been simplified significantly. Linux is not Unix, but it is a Unix-like operating system. Linux system is derived from Unix. Similar commands you will learn in this tutoral may be applicable to other unix systems.
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
-------

Linux is an operation system like Microsoft Windows or OSX. An operating system is a collection of software that manages system resources (such as memory, cpu(s), display and etc.) and provide applications a simpler interface to the system
hardware.

<img src="images/operating_system.png">
<img src="images/unix.png">
<img src="images/unix_like.png">
<img src="images/unix_family.png">
<img src="images/terminal.png">

## Before you start
-------
### Convention: 

***$*** denotes a prompt for the command line. It is not meant to be typed while running the commands in this tutorial. All Bash commands will appear in a gray box that you will run in your terminal. In this example below, you will only write "your\_command" in the terminal an press enter. It will then execute the command in your terminal. (Please Don't run "your\_command", since there is no such command :) )

	$ your_command

Let's start running some commands. First command will be "ssh" to connect to the UMASS cluster.

	$ ssh username@ghpcc06.umassrc.org
	

Let’s verify that we are at the right place.

	$ hostname
	
You need to see "ghpcc06" in your terminal. ghpcc06 is our "head node". You will learn what head or child nodes are in these tutorials later.

***Warning: Linux is case sensitive. If the command is "ls", you cannot run this command with LS.***

***Warning: There is no any spaces in the name of the commands in linux. However, if you need to use some parameters while running a command, you will use space in between a command and its parameters. For example; if you want to get a detailed list you will use ls -l and there is a space between ls and -l.***


## Getting started

To print something on the screen. We use echo command. 
	
	$ echo Hello World
	
If you just run 

	$ Hellow World
	-bash: Hello: command not found
	
You will get this error "command not found". Because there is no any command called "Hello".


## Basic Commands
<img src="images/dir_structure.png">



