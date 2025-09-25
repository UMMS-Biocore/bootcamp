# Preliminary Steps
1. [Step 1: Register for access to the HPC](#step-1-register-for-access-to-the-hpc-high-performance-computing-cluster)
2. [Step 2: Setup SSH Connection for your Computer](#step-2-setup-ssh-connection-for-your-computer)
    1. [If you're using Windows](#a-if-youre-using-windows)
    2. [If you're using Linux or MacOS](#b-if-youre-using-linux-or-macos)
3. [Step 3: Setup Connection for Via Foundry](#step-3-setup-connection-for-via-foundry)
3. [Step 4: Project Space Requirements](#step-4-project-space-requirements)

## Step 1: Register for access to the HPC (High Performance Computing Cluster)

The registration form can be found at <a href="https://hpc.umassmed.edu/wiki/index.php?title=Welcome_to_the_Scientific_Computing_for_Innovation_Cluster#Requesting_an_Account" target="_blank">HPC website</a> (If you're getting "**Operation timed out**" errors, try installing VPN software (eg. Pulse Secure) to access UMass Medical School network.
You can find the details at this <a href="https://umassmed.sharepoint.com/sites/information-technology/SitePages/VPN-Connect.aspx" target="_blank">UMass Medical School Link</a>.)

Once the HPC Admins receives your registration form, they will send an email to your PI requesting the PI’s permission to give you access. After it’s approved you will receive an email from the HPC Admins group with your HPC account user name.


## Step 2: Setup SSH Connection for your Computer

In order to use the pipelines in Via Foundry, you need access your cluster account. 

#### Troubleshooting: If you're getting "Operation timed out" errors, try installing VPN software (eg. Pulse Secure) to access UMass Medical School network.
You can find the details at this <a href="https://umassmed.sharepoint.com/sites/information-technology/SitePages/VPN-Connect.aspx" target="_blank">UMass Medical School Link</a>.

### A. If you're using Windows

In order to make an SSH connection to your account, you need to use program like PuTTY.

&nbsp;&nbsp;&nbsp;&nbsp;**A1.** Download and open PuTTY from <a href="https://www.putty.org/" target="_blank">their website</a>.

1. **Generate the Key:**
   - Run `puttygen` or `puttygen.exe`.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-1puttygensearch.png" width="50%">

   - Click the 'Generate' button.

   - Follow the instructions to move your mouse to generate randomness.
   - (Optional) Enter a passphrase in the 'Key passphrase' field and confirm it in the 'Confirm passphrase' field. This passphrase will be required each time you log into the HPC cluster.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-2generatekey.png" width="50%">

   - Click the 'Save public key' button and save the public key file (e.g., `public.key`). This file can be used later if needed.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-3savepublickey.png" width="50%">

   - Copy all the text in the box under "Public key for pasting into OpenSSH authorized_keys file" and paste it into the portal at [https://hpcportal.umassmed.edu/PublicKeys](https://hpcportal.umassmed.edu/PublicKeys). (Please use your UMASS email as username and email password for login.)
   
      Example public key:
      ```
      ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCHOGOPn5IaOL
      +yjA6KbIFVO5qoSq8rYWehXx9smUolajt5kGj71yEugchGs3BH
      ...
      Xv0QcmW9iFJTPxphFEH rsa-key-20240813
      ```

      - Please click "Create New" button
      - KeyName: Any name (e.g. My Laptop Key)
      - KeyValue: Public Key from your laptop.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/mac-terminal3.png" width="50%">

   - Click the 'Save private key' button and save the private key file (e.g., `private.ppk`). Remember the file name and location.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-3saveprivatekey.png" width="50%">

   - You can now close the `puttygen` application.

2. **Configure Putty for SSH Connection:**
   - Run `putty.exe` or the Putty application (not `puttygen`).

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-1puttysearch.png" width="50%">

   - If you don’t already have a session created for connecting to the SCI cluster, enter `hpc.umassmed.edu` as the hostname.
   - In the 'Saved Sessions' field, enter a name for the session (e.g., `HPC UMASS`) and click the 'Save' button to create a session for quick access.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-5savesession.png" width="50%">

      After Clicking on Save button:

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-5savesession2.png" width="50%">

   - Select your `HPC UMASS` session from the session list and click 'Load'.
   - In the Putty configuration, expand the 'Connection' category, then expand the 'SSH' category, and expand the 'Auth' category and then select 'Credentials' section.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-6credentials.png" width="50%">

   - In the 'Credentials' section, under 'Private key file for authentication', click 'Browse' and locate your private key file (e.g., `private.ppk`).

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/putty-6credentials2.png" width="50%">

   - After setting the location of the private key, go back to the 'Session' category, ensure that the `HPC UMASS` session is selected, and click 'Save'.
   - Now, when you launch the `HPC UMASS` session, Putty will use the configured public key. It will prompt you to enter the passphrase you set for the private key.

3. **Important Notes:**
   - **Username:** On the connection screen, enter your cluster username (e.g., `onur.yukselen-umw`).
   - **Passphrase:** When entering your passphrase, you won’t see any characters displayed as you type. Just type it and press Enter.
   - **Quick Access:** Since you configured your SSH connection session, each time you open Putty, you can simply select the `HPC UMASS` session, click 'Load', and then click 'Open' to connect.

### B. If you're using Linux or MacOS

If you are a MacOS or Linux user, you can make an SSH connection by using **Terminal**, a console program included with the operating system.

&nbsp;&nbsp;&nbsp;&nbsp;**B1.** Search for the **Terminal** program in your operating system and click to open it.
    
<img src="https://raw.githubusercontent.com/UMMS-Biocore/dolphinnext/master/docs/dolphinNext/dolphinnext_images/terminal.png" width="80%">
    
1. **Generate the SSH Key:**
   - From your command prompt, run the following command:
     ```
     ssh-keygen -t ecdsa -b 521
     ```
   - When prompted, accept the default file locations unless you already have keys with those names. Setting a password for the key file is strongly recommended.

2. **Locate Your Keys:**
   - You should end up with two keys in `~/.ssh` folder, use `cd ~/.ssh` command to enter that folder and `ls` command to see files:
      ```
      cd ~/.ssh
      ls
      ```

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/mac-terminal1.png" width="50%">

     - **Private key:** By default, this will be named `id_ecdsa`.
     - **Public key:** By default, this will be named `id_ecdsa.pub`.
     - Use `cat` command to print public key file content:

     ```
     cat id_ecdsa.pub
     ```

     <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/mac-terminal2.png" width="50%">

3. **Upload the Public Key:**
   - Open the following portal: [UMassMed HPC Public Keys](https://hpcportal.umassmed.edu/PublicKeys).
   - Copy the contents of your public key (`id_ecdsa.pub`) and paste it into the portal.

      - Please click "Create New" button
      - KeyName: Any name (e.g. My Laptop Key)
      - KeyValue: Public Key from your laptop.

         <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/mac-terminal3.png" width="50%">


4. **Connect with SSH command:**
   - Default approach to login your HPC account:
      ```
      ssh yourclusterusername@hpc.umassmed.edu
      ```
      
   - If you choose a name other than the default for your keys, you'll need to specify the location of your private key when authenticating via SSH.
       ```
       ssh -i /path/to/your/private/key yourclusterusername@hpc.umassmed.edu
       ```

## Step 3: Setup Connection for Via Foundry

1. Visit the Via Foundry website (https://viafoundry.umassmed.edu) and create SSH keys.
   - Click your profile avatar (top-right) → Credentials tab.
   - Click **Create credential**.
   - Set Credential type to **SSH** and enter any name.
   - Click **Generate keys for me**, then Save.
   - Next to your new key, click the **⋯ menu** → **View** → **Copy Button** to copy your public key. (You’ll paste this public key into the HPC website in the next step.)
  
2. Visit the HPC site and add SSH Keys to your account:
   - Add your public SSH key to the HPC site: https://hpcportal.umassmed.edu/PublicKeys/Create  (Please use your UMASS email as username and email password for login.)
   - Please click "Create New" button
   - KeyName: Any name (e.g. ViaFoundry Key)
   - KeyValue: Public Key from Foundry website.

      <img src="https://raw.githubusercontent.com/onuryukselen/bootcamp/master/images/mac-terminal3.png" width="50%">

3. At the Via Foundry website, Add a New Run environment for HPC
   - Click Profile Icon (at the top right) -> Click Run Environments Tab -> Click **Create Environment** Button
   - Choose the Profile Name "New UMASS SCI Cluster"
   - Enter your **HPC username** (you should receive an email from HPC Admins hpc@umassmed.edu. (example cluster username: john.thomas-umw))
   - Select your SSH Keys
   - Click "Test Connection" button to verify connection.
   - Click the "**Save**" button.


## Step 4: Project Space Requirements
Consult HPC-Admins (HPC-Admins@umassmed.edu) for your project space requirements. For example; typically 6 RNA-Seq libraries (5G to 10G each) require at least 500G of space to store the data and run the pipelines. Confirm you have the necessary space for your project. 
