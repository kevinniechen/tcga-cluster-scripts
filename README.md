### Genomic Data Commons
[Introduction to the GDC Portal](https://gdc.nci.nih.gov/access-data/gdc-data-portal)

[GDC Data Portal for TCGA Data](https://gdc-portal.nci.nih.gov/search/s)

[Guide to Downloading TCGA Data using the New GDC](https://gdc-docs.nci.nih.gov/Data_Transfer_Tool/Users_Guide/)

[Which GDC tools to use to download TCGA data](https://gdc.nci.nih.gov/access-data/data-access-processes-and-tools)

***

### Introduction to the ATRF Cluster
The NCI-ATRF cluster, which is a Moab HPC cluster, is a system of connected nodes, each representing a computer. When you first connect to the cluster *(see Sample Workflow: Connecting to the Cluster)*, you will be connected to the head node, which you can use for basic commands (ls, cd, mkdir), but you should not use it for any extensive scripts (otherwise it will slow down everyone else at NCI connected to the server on the head node).

Because of this, you should request a node to use in interactive mode `qsub -I` as soon as you connect, which essentially is allowing you to use a computer in the cluster that no one else is using. You can run small programs in this mode without doing anything else (./my-small-program), but this prevents you from being able to turn off your computer or closing the terminal without stopping whatever programs you ran in this node.

To allow programs to run 24/7 on the server regardless of whether your computer is on, you will need to submit a ***job*** to the server, which you can do with or without having started an interactive `qsub -I` node. To submit a job, you will first need to create a ***job file*** (call it myjobfile.pbs) in this format:

```bash
#!/bin/bash
#PBS -l nodes=1:ppn=1

cd $PBS_O_WORKDIR
echo 'Hello World' > output.txt
```
The first line is just a header to tell the interpreter to run bash. The second line starting with #PBS, is an option saying to run the job with 1 node and 1 processor. PBS options allow you to specify how you want the job to run (maximum time before it closes, how many processors to use, whether to email you if the job has finished, etc).

To submit the job file (*myjobfile.pbs*), enter in the command:
```bash
qsub myjobfile.pbs
```
The cluster will then begin running the commands specified in the job file. The Hello World job above will run and finish extremely quickly, but if you are running a file that takes a longer amount of time, you will probably want to check the status of your jobs at some point. To check the status of your running jobs, you can enter in the command:
```bash
pbsn -u username
```
Here is a real-world example of a job file that runs a bam-indexing script (not pictured here) and forwards the error and logs to specified directories:

```bash
#!/bin/bash
#PBS -N bam-index
#PBS -q medium   
#PBS -l nodes=1
#PBS -l cput=24:00:00        
#PBS -m bea
#PBS -M chenk8@mail.nih.gov

cd $PBS_O_WORKDIR
./bam-index.sh > ~/err/$PBS_JOBID > ~/logs/$PBS_JOBID
exit 0
```

[Further explanation on the various job script options can be found here.](https://hpcc.usc.edu/support/documentation/running-a-job-on-the-hpcc-cluster-using-pbs/)

***

### Sample workflow
#### Login to the Moab HPC cluster
Open your terminal and log in with your NIH username and password:
```bash
ssh USERNAME@moab.ncifcrf.gov
```

Then, open up an interactive job (so that the commands/scripts that you run
won't slow down the head node of the server):
```bash
qsub -I
```
Enter in `exit` if you ever want
to leave it.

Change to the lab directory:
```bash
cd /ifs/projects/GRCBL-NGS/slowtemp
```

#### Start a project
Create a project directory:
```bash
mkdir project1/
cd project1/
```

Create the cghub key:
```bash
nano cghub.key
```
Then copy and paste your TCGA key (in your email) into the folder, and enter Ctrl+X and then Y
to save.

Download basic TCGA and BAM scripts and set them as executable:
```bash
git clone https://github.com/kevchn/tcga-cluster-scripts scripts/
chmod +x scripts/*.sh
```

#### Getting your TCGA authentication token
You will need an authentication token to download controlled data from TCGA (most raw-sequencing files). Minimize the terminal (but do not close it).

To generate a token, first log in to the [GDC Data Portal](https://gdc-portal.nci.nih.gov/) by clicking the Login button in the top right corner of the page.

After logging in, clicking the username will open a drop-down menu. Select Download Token from the menu to generate an authentication token. Copy the contents of this token and save it `nano token.txt` in the project1/ directory.

#### Finding the UUID of a single TCGA sample
Head to the [Genomic Data Commons Data Portal](https://gdc-portal.nci.nih.gov/search/s)

Click on any specific file to get more information, including its UUID (for example: 22a29915-6712-4f7a-8dba-985ae9a1f005)

#### Downloading a single TCGA file onto the server
If you have the UUID of a single TCGA file, use the new [gdc-data-transfer-tool](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool):
```bash
module add gdc-client
gdc-client download 22a29915-6712-4f7a-8dba-985ae9a1f005 -t token.txt
```
Please note that the NCI-ATRF cluster may not have the gdc-client installed, or it may be under a different name, in which case the above commands will not work. Enter in `module list` to see if there are any modules called gdc, genomic-data-commons, or similar, and import that module instead. If no modules exist, contact the NCI helpdesk to request that the program be installed.

[You can also download the gdc-client locally on your computer](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool).

#### Download a list of desired TCGA files
Minimize the terminal (but do not close it), and head to the [Genomic Data Commons Data Portal](https://gdc-portal.nci.nih.gov/search/s)

Use the sidebar to filter out the desired samples. For example, if you want all
GBM RNA-seq samples, under the Cases tab, select Primary Site>Brain and Project>TCGA-GBM, and then under the Files tab, select Data Category>Raw Sequencing Data and Experimental Strategy>RNA-seq. The end-result should be 174 files.

Click the Download Manifest button to download an xml file containing the IDs of all these files (so that you can download them onto the cluster). [More information regarding how to download individual and multiple IDs from TCGA can be found here.](https://gdc-docs.nci.nih.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/)

Open the xml file, copy the contents, and then open up terminal again.

In your project1/ folder, `nano manifest.xml` and paste in the XML information before saving (Ctrl+X, then Y).

#### Bulk downloading TCGA files onto the server
Once you have the manifest.xml file, you can download the files using the new [gdc-data-transfer-tool](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool):
```bash
module add gdc-client
gdc-client download -m manifest.xml -t token.txt
```

You will most likely want to do this using a job script, so that the cluster will download the files 24/7 until the files are all downloaded.
***
### Running this repository's jobfiles and scripts
The following usage commands assume that you are in the parent directory of the scripts/ folder. If you are in a different directory, for example, projects1/data, you will have to provide either a different relative path for usage: ```qsub ../scripts/index-bams.pbs``` or an absolute path: ```qsub /ifs/projects/GRCBL-NGS/slowtemp/project1/scripts/index-bams.pbs```.

#### scripts/index-bams.pbs
Indexes (and possibly sorts) all .bam files in the current directory

Creates .bai files for all .bam files in the current directory. Doesn't change the original .bam files. Need to ensure that the bam files are already sorted. TCGA files are usually sorted. To turn on sorting befor indexing, ```nano scripts/index-bams.pbs``` and remove the '#' from line 15 and 16 and comment out line 17 like below:
```bash
samtools sort "${bamfile}" "${bamfile}.sorted"
samtools index ${bamfile}.sorted > err/$PBS_JOBID > logs/$PBS_JOBID
# samtools index $bamfile > err/$PBS_JOBID > logs/$PBS_JOBID
```

Usage:
* ```qsub scripts/index-bams.pbs```

Input:
* sample1.bam
* sample2.bam
* sample3.bam

Returns:
* sample1.bai
* sample2.bai
* sample3.bai

#### scripts/separate-dif-len-bams.pbs
Moves .bam files to new directories based on read-length and creates a summary of .bam file read-lengths

Checks the read-length of the .bam files by looking at the header of the bam files and then moves the files into a new directory based on that read-length. Removes bam-files with multiple read-lengths. Used for programs that can only run on multiple .bam files of the same read-length. Summary file, bam-information.txt, gives read-count and read-length of each bam-file.

Usage:
* ```qsub scripts/separate-dif-len-bams.pbs```

Input:
* sample1.bam
* sample2.bam
* sample3.bam

Returns:
* bam-information.txt
* bams-50/sample1.bam
* bams-50/sample2.bam
* bams-75/sample3.bam
