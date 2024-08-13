## Step 1: Register for access to the HPC (High Performance Computing Cluster)

The registration form can be found at <a href="https://hpc.umassmed.edu/wiki/index.php?title=Welcome_to_the_Scientific_Computing_for_Innovation_Cluster#Requesting_an_Account" target="_blank">HPC website</a>. Once the HPC Admins receives your registration form, they will send an email to your PI requesting the PI’s permission to give you access. After it’s approved you will receive an email from the HPC Admins group with your HPC account user name.


## Step 2: Connecting to your account 

In order to use the pipelines in ViaFoundry/DolphinNext, you need access your cluster account. 

#### Troubleshooting: If you're getting "Operation timed out" errors, try installing VPN software (eg. Pulse Secure) to access UMass Medical School network.
You can find the details at this <a href="https://umassmed.sharepoint.com/sites/information-technology/SitePages/VPN-Connect.aspx" target="_blank">UMass Medical School Link</a>.

### A. If you're using Windows

In order to make an SSH connection to your account, you need to use program like PuTTY.

&nbsp;&nbsp;&nbsp;&nbsp;**A1.** Download and open PuTTY from <a href="https://www.putty.org/" target="_blank">their website</a>.

1. **Generate the Key:**
   - Run `puttygen` or `puttygen.exe`.
   - Click the 'Generate' button.
   - Follow the instructions to move your mouse to generate randomness.
   - (Optional) Enter a passphrase in the 'Key passphrase' field and confirm it in the 'Confirm passphrase' field. This passphrase will be required each time you log into the HPC cluster.
   - Click the 'Save public key' button and save the public key file (e.g., `public.key`). This file can be used later if needed.
   - Copy all the text in the box under "Public key for pasting into OpenSSH authorized_keys file" and paste it into the portal at [https://hpcportal.umassmed.edu/PublicKeys](https://hpcportal.umassmed.edu/PublicKeys).
   
   Example public key:
   ```
   ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCHOGOPn5IaOL
   +yjA6KbIFVO5qoSq8rYWehXx9smUolajt5kGj71yEugchGs3BH
   dvE3zptIGbLt4uXRyJxb4JtgBOqnYq43o3AeFGhqSfcamWid/d
   1IbXr7Ii6gYmGKJwquIGU9d29IWHvLaFICnxZKFXxtsJRxZcc0
   XLN1eKRxz/nj3jIMUIG1iFYelyrk6I4nZ0zcBYGFTt76xln1Yb
   QCehM0fOFhMw2xyuxT8tGfixHSc+b0Lcie7UijYPAB+G9mMKzR
   bQaBFzbeX9ecyU1dUTM1WrgbNKavGXv0QcmW9iFJTPxphFEH rsa-key-20240813
   ```

   - Click the 'Save private key' button and save the private key file (e.g., `private.ppk`). Remember the file name and location.
   - You can now close the `puttygen` application.

2. **Configure Putty for SSH Connection:**
   - Run `putty.exe` or the Putty application (not `puttygen`).
   - If you don’t already have a session created for connecting to the SCI cluster, enter `hpc.umassmed.edu` as the hostname.
   - In the 'Saved Sessions' field, enter a name for the session (e.g., `HPC UMASS`) and click the 'Save' button to create a session for quick access.
   - Select your `HPC UMASS` session from the session list and click 'Load'.
   - In the Putty configuration, expand the 'Connection' category, then expand the 'SSH' category, and select the 'Auth' category.
   - In the 'Credentials' section, under 'Private key file for authentication', click 'Browse' and locate your private key file (e.g., `private.ppk`).
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
    
- From the command prompt, please run 'ssh-keygen -t ecdsa -b 521' and accept the default locations for the file locations unless you already have keys with those names. Setting a password for the key file is strongly encouraged.
- You should end up with two keys, the private one (defaults to id_ecdsa) and the public one (defaults to id_ecdsa.pub). Please paste the contents of the public key into the portal, https://hpcportal.umassmed.edu/PublicKeys
- If you choose a name other than the default for your keys, you will need to tell ssh where to find the private key to authenticate. You can do this by specifying "-i /path/to/your/private/key", or by creating a ~/.ssh/config file and adding a Host entry, e.g.:
	
```
ssh -i ~/.ssh/id_ecdsa yourclusterusername@hpc.umassmed.edu
```

## Step 3: Setup Connection for Dolphinnext/Viafoundry

1. Visit the Via Foundry website (https://viafoundry.umassmed.edu) and create SSH keys.
   - Click Profile Icon (at the top right) -> Click SSH Keys Tab -> Click Add SSH Key Button
   - Enter any name for your SSH key.
   - Click the "Create new keys" checkbox and click the "Generate Keys" button.
   - Copy your public SSH key and click the "Submit" button.
2. Visit the HPC site and add SSH Keys to your account:
   - Add your public SSH key to the HPC site: https://hpcportal.umassmed.edu/PublicKeys/Create
   - Note: Please use your UMASS email as username and email password for login.
3. At the Via Foundry website, Add a New Run environment for HPC
   - Click Profile Icon (at the top right) -> Click Run Environments Tab -> Click Add Environment Button
   - Choose the Profile Name "New UMASS SCI Cluster"
   - Enter your HPC username
   - Select your SSH Keys
   - Click "Check Connection" button to verify connection.
   - Click the "**Save Changes**" button.


## Step 4: Project Space Requirements
Consult HPC-Admins (HPC-Admins@umassmed.edu) for your project space requirements. For example; typically 6 RNA-Seq libraries (5G to 10G each) require at least 500G of space to store the data and run the pipelines. Confirm you have the necessary space for your project. 
