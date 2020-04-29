*****************
Quick Start Guide
*****************

Signing Up
==========

This guide will walk you through how to start using DolphinNext pipelines. First off, you need to enter DolphinNext web page: https://dolphinnext.umassmed.edu/ and click **Sign Up** or **Sign in with Google** buttons. You will be asked to enter your institute information. An email will be sent to you upon verification of your information. 

.. image:: dolphinnext_images/sign_in.png
	:align: center
	:width: 99%

Creating Profile
================

Once you enter DolphinNext platform, **Profile Wizard** will open as shown at below. 

.. image:: dolphinnext_images/profilewizard_first.png
	:align: center
	:width: 99%

.. note::
    If you need to re-open wizard window, you can simply click wizard button at the top right.

    .. image:: dolphinnext_images/profilewizard_button.png
	   :align: center
	   :width: 45%

.. note:: If you have any issues/questions about creating profiles please contact us on: support@dolphinnext.com

1. First, choose your profile type:

    A. If you have an access to High Performance Computing (HPC) environments, or personal workstations, then please choose **Host**. 
    B. If you have an Amazon Web Services (AWS) account or planning to create one, then please choose **Amazon** and follow our `Amazon Guide <https://dolphinnext.readthedocs.io/en/latest/dolphinNext/profile.html#b-defining-amazon-profile>`_  to create your run environment.
    C. If you have an Google Cloud account or planning to create one, then please choose **Google** and follow our `Google Guide <https://dolphinnext.readthedocs.io/en/latest/dolphinNext/profile.html#c-defining-google-profile>`_ to create your run environment.
    D. If you choose **MGHPCC cluster** option, you can upload your files to our MGHPCC cluster to process your data and download your results from report section. However, you will not have direct access to our cluster.
    
.. note::  If you choose MGHPCC option, you can skip the rest of this guideline and go to `Running Pipelines Section <quick.html#running-pipelines>`_.
        
2. Second, add public SSH Key into your host machine.

    -  Please confirm our `Terms and Conditions <https://dolphinnext.umassmed.edu/php/terms.php>`_ & `Privacy Policy <https://dolphinnext.umassmed.edu/php/privacy.php>`_ by clicking checkbox.
    -  Here, public key is securely generated for your account and required to be added into ``~/.ssh/authorized_keys`` in the host by user. Please check our `Adding Public SSH Key Section <public_ssh_key.html>`_ for help. 
    -  After adding public key, please click Validate SSH Keys button to finalize this section.
    
.. important::  **Username/Hostname:** You should enter your username and hostname of the host which you would like to connect (yourusername@yourhostname). For instance, for us2r@ghpcc06.umassrc.org::
    
        *  Username: yourusername (eg. us2r)
        *  Hostname: yourhostname (eg. ghpcc06.umassrc.org)
    

3. Third, install/validate software dependencies into the host machine.

In order to execute our pipelines, Nextflow should be installed into your host environment. Besides, most of our pipelines isolates their dependencies within their Docker or Singularity containers, therefore please install these softwares into your machine by following guidelines. If your platform doesn't support the installation of Docker, you can still use our pipelines by just using Singularity.
    
    -  Installing `Nextflow <faq.html#how-can-i-install-nextflow>`_ 
    -  Installing `Singularity (Version 3) <faq.html#how-can-i-install-singularity>`_ 
    -  Installing `Docker <faq.html#how-can-i-install-docker>`_ 

    
    Sofware Dependencies Section:

    -  **JAVA Command (optional):** If JAVA is not added to $PATH environment, you can run command (eg. ``module load java/8.0``) to manipulate your $PATH environment and gain access to JAVA.
    
    -  **Nextflow Path or Command (optional):** If nextflow is not added to $PATH environment, you can either enter the path of the nextflow (eg. ``/project/bin``), or run command (eg. ``module load nextflow``) to manipulate your $PATH environment and gain access to new softwares. 
    
    -  **Docker/Singularity Command (optional):** You can run command to manipulate your $PATH environment in order to gain access to new softwares. (eg. ``module load docker/1.0.0`` or ``module load singularity/3.0.0``) 
    
4. General run settings could be set in following **Run settings** section:

    -  **Executor of Nextflow:** Nextflow itself is initiated with this method and it will be only used for running nextflow itself.
    
    -  **Executor of Nextflow Jobs:** This setting will be used as default setting for submitted jobs by Nextflow.
    
    -  **Download Directory:** Used to download shared pipeline files such as genome indexes. If your platform has already such path, please enter that location. Otherwise you can set any path that you have permission to write. e.g. ``/share/dolphinnext/downloads``

Once you complete these steps, you're now able to start using publicly available pipelines.


Running Pipelines
=================

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/gaq_LwewFPA" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
    </br>


1. The easiest way to run pipeline is using main page by clicking the **Biocore DolphinNext** button at the top left of the screen. Now, you can investigate publicly available pipelines as shown at below and select the pipeline you want run by clicking **Learn More** button.

    .. image:: dolphinnext_images/main_page.png
	   :align: center


2. Once pipeline is loaded, you will notice "Run" button at the right top of the page.


    .. image:: dolphinnext_images/project_runbutton.png
	   :align: center
	   :width: 35%


3. This button opens new window where you can create new project by clicking "Create a Project" button. After entering and saving the name of the project, it will be added to your project list. Now you can select your project by clicking on the project as shown in the figure below.

    .. image:: dolphinnext_images/project_pipeselect.png
	   :align: center

4. Now, you may proceed with entering run name which will be added to your run list of the project. Clicking "Save run" will redirect to "run page".

5. Initially, in the header of the run page, orange ``Waiting`` button will be shown. In order to initiate run, following data need to be entered:

    .. image:: dolphinnext_images/run_header_waiting.png
	   :align: center

    A. **Work Directory:**  Full path of the directory, where nextflow runs will be executed.
    
        .. image:: dolphinnext_images/run_params_work.png
	   :align: center
	   :width: 99%
    
    B. **Run Environment:** Profile that is created in the `profile <profile.html>`_  page. If `Amazon profile <profile.html#b-defining-amazon-profile>`_ or `Google profile <profile.html#c-defining-google-profile>`_ is selected, then status of the profile should to be at the stage of **running**.
    
        .. image:: dolphinnext_images/run_params_env.png
	   :align: center
	   :width: 99%
    
    C. **Inputs:** Value and path of the files need to be entered. For detailed information please check `adding files section. <quick.html#adding-files>`_ 

        .. image:: dolphinnext_images/run_params_inputs.png
	   :align: center
	   :width: 50%


6. Once all requirements are satisfied, ``Waiting`` button will turn in to green ``ready to run`` button as shown below. You can initiate your run by clicking ``ready to run`` button. Please go through `run page <run.html>`_ for detailed explanation about each module is used.
    
    .. image:: dolphinnext_images/run_header_ready.png
	   :align: center



Adding Files
============

Remote Files
------------
You can reach your remote files by entering:

    - Full path of a directory: eg. ``/share/data/umw_biocore/genome_data/mousetest/mm10/gz``
    - Web link: eg. ``https://galaxyweb.umassmed.edu/pub/dnext_data/test/reads``
    - Amazon (S3) Bucket: eg. ``s3://biocore/fastq``
    - Google (GS) Bucket: eg. ``gs://biocore/fastq``

Geo Files
---------

If you want to download and use NCBI (GEO data) in the pipeline, you can simply use **GEO Files** tab. Here are the few examples for GEO ID: ``GSM1331276``, ``GSE55190``, ``SRR10095965``
    
Upload Files
------------
If you need to upload your local files and transfer into **Target Directory in the Host**, you can use **Upload Files** tab.

For detailed information about adding files, please check our tutorial video:
        
    .. raw:: html

        <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe src="https://www.youtube.com/embed/3QaAqdyB11w" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
        </div>
        </br>
 
 
 





