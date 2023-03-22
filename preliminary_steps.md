# Introduction


  1. Connecting cluster (Authentication Section):

  2. Setup Connection for Dolphinnext/Viafoundry
1. Visit the Via Foundry website (https://viafoundry.umassmed.edu) and create SSH keys.
   a. Click Profile Icon (at the top right) -> Click SSH Keys Tab -> Click Add SSH Key Button
   b. Enter any name for your SSH key.
  c. Click the "Create new keys" checkbox and click the "Generate Keys" button.
   d. Copy your public SSH key and click the "Submit" button.
2. Visit the HPC site and add SSH Keys to your account:
   a. Add your public SSH key to the HPC site: https://hpcportal.umassmed.edu/PublicKeys/Create
   Note: Please use your UMASS email as username and email password for login.
3. At the Via Foundry website, Add a New Run environment for HPC
   a. Click Profile Icon (at the top right) -> Click Run Environments Tab -> Click Add Environment Button
   b. Choose the Profile Name "New UMASS SCI Cluster"
   c. Enter your HPC username
   d. Select your SSH Keys
   e. Click the "save changes" button.
If you've any questions about these steps, please let me know.
Best,
Onur

# Introduction

Please follow these step to begin using Dolphin services:

  * [Step 1: Register for access to the HPC](#step-1-register-for-access-to-the-hpc-high-performance-computing-cluster)
  * [Step 2: Send an email to the HPCC-Admins](#step-2-send-an-email-to-the-hpcc-admins)
  * [Step 3: Connecting to your account](#step-3-connecting-to-your-account)
    * [A. If you're using Windows](#a-if-youre-using-windows)
    * [B. If you're using Linux or MacOS](#b-if-youre-using-linux-or-macos)
  * [Step 4: Run a script for authorization](#step-4-run-a-script-for-authorization)
  * [Step 5: Project Space Requirements](#step-5-project-space-requirements)
  

## Step 1: Register for access to the HPC (High Performance Computing Cluster)

The registration form can be found at <a href="https://hpc.umassmed.edu/wiki/index.php?title=Welcome_to_the_Scientific_Computing_for_Innovation_Cluster#Requesting_an_Account" target="_blank">HPC website</a>. Once the HPC Admins receives your registration form, they will send an email to your PI requesting the PI’s permission to give you access. After it’s approved you will receive an email from the HPC Admins group with your HPC account user name.


## Step 2: Connecting to your account 

In order to use the pipelines in ViaFoundry/DolphinNext, you need access your cluster account. 

#### Troubleshooting: If you're getting "Operation timed out" errors, try installing VPN software (eg. Pulse Secure) to access UMass Medical School network.
You can find the details at this <a href="https://umassmed.sharepoint.com/sites/information-technology/SitePages/VPN-Connect.aspx" target="_blank">UMass Medical School Link</a>.

### A. If you're using Windows

In order to make an SSH connection to your account, you need to use program like PuTTY.

&nbsp;&nbsp;&nbsp;&nbsp;**A1.** Download and open PuTTY from <a href="https://www.putty.org/" target="_blank">their website</a>.

- Run 'puttygen.exe'
- Click on the 'Generate' button.
- Follow instructions to move your mouse.
- Enter a passphrase that you will remember in the 'Key passphrase' section. Enter the same passphrase in the 'Confirm passphrase' section.
- Copy all the text in the box under "Public key for pasting into OpenSSH authorized_keys file and paste into the portal, https://hpcportal.umassmed.edu/PublicKeys
- Click the 'Save private key' button and save the private key file. Remember the name and location of where you save the private key.
- Run 'putty.exe'
- If you don't already have a session created for connecting to the SCI cluster, enter 'hpc.umassmed.edu' as the hostname and saved session name and click the 'Save' button to create one.
- Select your 'hpc.umassmed.edu' session from the session list and select 'Load'.
- In the putty configuration for your saved SCI connection expand the 'Connection' category, then the 'SSH' category and then select the 'Auth' category.
- Under 'Private key file for authentication' you can enter the location and name of your private key, or select the 'Browse' button to locate the private key.
- Once you've set the location of the private key, go back to the 'Session' category, make sure 'hpc.umassmed.edu' is still selected, and press the 'Save' button.
- Now when you launch the 'hpc.umassmed.edu' session putty knows the public key to try. It will prompt you for the passphrase you set for the private key.

### B. If you're using Linux or MacOS

If you are a MacOS or Linux user, you can make an SSH connection by using **Terminal**, a console program included with the operating system.

&nbsp;&nbsp;&nbsp;&nbsp;**B1.** Search for the **Terminal** program in your operating system and click to open it.
    
<img src="https://raw.githubusercontent.com/UMMS-Biocore/dolphinnext/master/docs/dolphinNext/dolphinnext_images/terminal.png" width="80%">
    
- From the command prompt, please run 'ssh-keygen -t ecdsa -b 521' and accept the default locations for the file locations unless you already have keys with those names. Setting a password for the key file is strongly encouraged.
- You should end up with two keys, the private one (defaults to id_ecdsa) and the public one (defaults to id_ecdsa.pub). Please paste the contents of the public key into the portal, https://hpcportal.umassmed.edu/PublicKeys
- If you choose a name other than the default for your keys, you will need to tell ssh where to find the private key to authenticate. You can do this by specifying "-i /path/to/your/private/key", or by creating a ~/.ssh/config file and adding a Host entry, e.g.:
	
```
ssh -i ~/.ssh/id_ecdsa yourclusterusername@hpc.umassmed.edu
```

## Step 3: Setup Connection for Dolphinnext/Viafoundry

1. Visit the Via Foundry website (https://viafoundry.umassmed.edu) and create SSH keys.
   a. Click Profile Icon (at the top right) -> Click SSH Keys Tab -> Click Add SSH Key Button
   b. Enter any name for your SSH key.
  c. Click the "Create new keys" checkbox and click the "Generate Keys" button.
   d. Copy your public SSH key and click the "Submit" button.
2. Visit the HPC site and add SSH Keys to your account:
   a. Add your public SSH key to the HPC site: https://hpcportal.umassmed.edu/PublicKeys/Create
   Note: Please use your UMASS email as username and email password for login.
3. At the Via Foundry website, Add a New Run environment for HPC
   a. Click Profile Icon (at the top right) -> Click Run Environments Tab -> Click Add Environment Button
   b. Choose the Profile Name "New UMASS SCI Cluster"
   c. Enter your HPC username
   d. Select your SSH Keys
   e. Click "Check Connection" button to verify connection.
   e. Click the "**save changes**" button.


## Step 4: Project Space Requirements
Consult HPC-Admins (HPC-Admins@umassmed.edu) for your project space requirements. For example; typically 6 RNA-Seq libraries (5G to 10G each) require at least 500G of space to store the data and run the pipelines. Confirm you have the necessary space for your project. 
