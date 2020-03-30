# Introduction

Please follow these step to begin using Dolphin services:

  * [Before you start](#step-1-register-for-access-to-the-hpcc)
  * [Getting Started](#getting-started)
  * [Exercise 1: Creating processes](#exercise-1---creating-processes)
    * [FastQC process](#1-fastqc-process)
    * [Hisat2 process](#2-hisat2-process)
    * [RSeQC process](#3-rseqc-process)
  * [Exercise 2: Building a pipeline](#exercise-2---building-a-pipeline)
  * [Exercise 3: Running a pipeline](#exercise-3---executing-a-pipeline)
  * [Exercise 4 (optional): Supporting single and paired end reads](#exercise-4-optional---supporting-both-single-and-paired-end-reads)

## Step 1: Register for access to the HPCC (High Performance Computing Cluster)

The registration form can be found at <a href="https://www.umassrc.org/hpc/" target="_blank">MGHPC website</a>. Once the HPCC Admins receives your registration form, they will send an email to your PI requesting the PI’s permission to give you access. After it’s approved you will receive an email from the HPCC Admins group with your HPCC account user name.

## Step 2: Send an email to the HPCC-Admins 

Send an email to the HPCC-Admins (HPCC-Admins@umassmed.edu) to set your account for Dolphin.

## Step 3: Connecting to your account 

In order to use the pipelines in DolphinNext, you need to run a script in the cluster. This is a one time script that will allow DolphinNext to submit future jobs to the cluster on your behalf. To connect cluster use your HPCC user ID and password sent by HPCC Admins.

### A. If you're using Windows

In order to make an SSH connection to your account, you need to use program like PuTTY.

&nbsp;&nbsp;&nbsp;&nbsp;**A1.** Download and open PuTTY from <a href="https://www.putty.org/" target="_blank">their website</a>.

&nbsp;&nbsp;&nbsp;&nbsp;**A2.** Use the following info to configure your connection and click **Open** to start the SSH session.
    
    - Host Name: ghpcc06.umassrc.org
    - Port: 22 
    - Connection Type: SSH 
        
<img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty_ghpcc.png" width="50%">

&nbsp;&nbsp;&nbsp;&nbsp;**A2.** If this is the first time that you have used PuTTY to log in to your account with SSH, a warning similar to the following one displays. If you are sure that you have entered the correct information, click Yes. 

<img src="https://raw.githubusercontent.com/UMMS-Biocore/dolphinnext/master/docs/dolphinNext/dolphinnext_images/putty_warning.png" width="40%">

&nbsp;&nbsp;&nbsp;&nbsp;**A4.** After you accept the warning, the terminal prompts you for your username and password. Please enter these values and press enter.

```
Caution:
* The password is not echoed back to the screen as it is entered.
* If you need to copy and paste your password, you can right-click (or use middle mouse button) to paste your password.
```

&nbsp;&nbsp;&nbsp;&nbsp;**A5.** If this is the first time that you login to your account, you might need to reset your password. Please enter the new password and press Enter. Note that passwords are not echoed to the screen. After changing your password, you will be logged out and you need to reconnect to your host machine by re-opening Putty window.

&nbsp;&nbsp;&nbsp;&nbsp;**A6.** If you have entered the correct password, the prompt responds with a shell prompt::

	[yourusername@ghpcc06 ~]#

### B. If you're using Linux or MacOS

If you are a MacOS or Linux user, you can make an SSH connection by using **Terminal**, a console program included with the operating system.

&nbsp;&nbsp;&nbsp;&nbsp;**B1.** Search for the **Terminal** program in your operating system and click to open it.
    
<img src="https://raw.githubusercontent.com/UMMS-Biocore/dolphinnext/master/docs/dolphinNext/dolphinnext_images/terminal.png" width="80%">
    
&nbsp;&nbsp;&nbsp;&nbsp;**B2.** First, you should type the ``ssh`` command in the console. Then enter your username and hostname and add ``@`` sign in between and press enter::
    
        ssh yourusername@ghpcc06.umassrc.org
        
<img src="https://raw.githubusercontent.com/UMMS-Biocore/dolphinnext/master/docs/dolphinNext/dolphinnext_images/terminal_ssh.png" width="80%">

&nbsp;&nbsp;&nbsp;&nbsp;**B3.** The terminal prompts you for your password. Please enter your password and press enter.

```
Caution:
* The password is not echoed back to the screen as it is entered.
* You can still **copy and paste your password**, but it won't appear in your screen.
```

&nbsp;&nbsp;&nbsp;&nbsp;**B4.** If this is the first time that you login to your host machine, you might need to reset your password. Please enter the new password and press Enter. Note that passwords are not echoed to the screen. After changing your password, you will be logged out and you need to reconnect to your account.

&nbsp;&nbsp;&nbsp;&nbsp;**B5.** If you have entered the correct password, the prompt responds with a shell prompt::

        [yourusername@ghpcc06.umassrc.org ~]#

## Step 4: Run a script for authorization

In order to use the pipelines in DolphinNext, you need to run a script in the cluster. This is a one time script that will allow DolphinNext to submit future jobs to the cluster on your behalf. 

#### 1. Please enter following command to your terminal and press enter.

```
/project/umw_biocore/bin/addKey.bash
```

#### 2. Command will print following text: 
```
[yourusername@ghpcc06 ~]$ /project/umw_biocore/bin/addKey.bash
ummsres15 key added
ummsres03 key added
```    

#### 3. Send an email to the Biocore

Please send “the output” of this script and your cluster username to Biocore (biocore@umassmed.edu). We will make sure that you successfully added the keys to your cluster system.  

## Step 5: Project Space Requirements
Consult HPCC-Admins for your project space requirements. For example; typically 6 RNA-Seq libraries (5G to 10G each) require at least 500G of space to store the data and run the pipelines. Confirm you have the necessary space for your project. 
