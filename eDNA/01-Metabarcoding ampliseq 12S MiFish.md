# Metabarcoding workflow for 12S amplicon sequencing with MiFish primers

*page details in progress.* 

The 12S rRNA gene region of the mitogenome is ~950 bp. There are two popular primer sets to amplify two different regions of 12S: Riaz and MiFish. The following workflow includes script specific to the **MiFish Universal (U) or Elasmobranch (E)** primer set.

![](https://th.bing.com/th/id/OIP.EbXPYETYLBPEymNEIVEGLQHaCc?rs=1&pid=ImgDetMain)

Citation: [Miya et al. 2015](https://royalsocietypublishing.org/doi/full/10.1098/rsos.150088)

The MiFish-U and MiFish-E primers are two variants of universal PCR primers developed for metabarcoding environmental DNA (eDNA) from fishes. Here are the key differences between them:  
- *Target species*: MiFish-U (Universal) is designed to amplify DNA from a wide range of bony fishes (Osteichthyes). MiFish-E (Elasmobranch) is specifically optimized for cartilaginous fishes like sharks and rays (Elasmobranchii).   
- *Primer sequences*: While both primer sets target the mitochondrial 12S rRNA gene, they have slightly different nucleotide sequences to accommodate the genetic variations between bony and cartilaginous fishes.  
- *Amplicon length*: MiFish-U typically produces amplicons around 170-180 base pairs long and MiFish-E amplicons are usually slightly shorter, around 160-170 base pairs. 

Scripts to run: 

1. 00-fastqc.sh   
2. 00-multiqc.sh  
3. 01a-metadata.R
4. 01b-ampliseq.sh
5. 02-taxonomicID.sh  

## Step 1: Conda environment: Fisheries eDNA 

GMGI Fisheries has a conda environment set-up with all the packages needed for this workflow. Code below was used to create this conda environment. **DO NOT REPEAT** every time user is running this workflow.

```
# Activate conda
source ~/../../work/gmgi/miniconda3/bin/activate

# Creating conda 
conda create --name fisheries_eDNA

# Installing packages needed for this workflow 
conda install -c bioconda fastqc 
conda install multiqc 
conda install bioconda::nextflow 
conda install conda-forge::singularity
conda install bioconda::blast
conda install nextflow
```

The conda environment is started within each slurm script, but to activate conda environment outside of the slurm script to update packages or check what is installed:

```
# Activate conda
source ~/../../work/gmgi/miniconda3/bin/activate

# Activate fisheries eDNA conda environment 
conda activate fisheries_eDNA

# List all available environments 
conda env list 

# List all packages installed in fisheries_eDNA
conda list 

# Update a package
conda update [package name]
``` 
 
## Step 2: Assess quality of raw data  

`00-fastqc.sh`: 

```
#!/bin/bash
#SBATCH --error=script_output/fastqc_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=script_output/fastqc_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=fastqc
#SBATCH --mem=3GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for your project

# Activate conda environment
source ~/../../work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

## LOAD MODULES
module load OpenJDK/19.0.1 ## dependency on NU Discovery cluster 
module load fastqc/0.11.9

## SET PATHS 
raw_path=""
out_dir=""

## CREATE SAMPLE LIST FOR SLURM ARRAY
### 1. Create list of all .gz files in raw data path
ls -d ${raw_path}/*.gz > ${raw_path}/rawdata

### 2. Create a list of filenames based on that list created in step 1
mapfile -t FILENAMES < ${raw_path}/rawdata

### 3. Create variable i that will assign each row of FILENAMES to a task ID
i=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

## RUN FASTQC PROGRAM 
fastqc ${i} --outdir ${out_dir}
```

To run:    
- Start slurm array (e.g., with 138 files) = `sbatch --array=0-136 01-fastqc.sh`.

Notes:  
- This is going to output *many* error and output files. After job completes, use `cat *output.* > ../fastqc_output.txt` to create one file with all the output and `cat *error.* > ../fastqc_error.txt` to create one file with all of the error message outputs. 
- Within the `out_dir` output folder, use `ls *html | wc` to count the number of html output files (1st/2nd column values). This should be equal to the --array range used and the number of raw data files. If not, the script missed some input files so address this before moving on.  


## Step 3: Visualize quality of raw data  

`00-multiqc.sh` 

```
#!/bin/bash
#SBATCH --error=output_messages/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=multiqc
#SBATCH --mem=8GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for your project
## 2. Optional: change file name (multiqc_raw.html) as desired

# Activate conda environment
source ~/../../work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

## SET PATHS 
## fastqc_output = output from 00-fastqc.sh; fastqc program
fastqc_output=""
multiqc_dir=""

## RUN MULTIQC 
multiqc --interactive ${fastqc_output} -o ${multiqc_dir} --filename multiqc_raw.html
```

To run:  
- `sbatch 00-multiqc.sh` 

Notes:  
- Depending on the number of files per project, multiqc can be quick to run without a slurm script. To do this, run each line separately in the command line after activating the conda environment.  

## Step 4: nf-core/ampliseq 

#### 12S primer sequences (required)

Below is what we used for 12S amplicon sequencing. Ampliseq will automatically calculate and include the reverse compliment sequence. 

MiFish-U 12S amplicon F: GTCGGTAAAACTCGTGCCAGC  
MiFish-U 12S amplicon R: CATAGTGGGGTATCTAATCCCAGTTTG       

MiFish-E 12S amplicon F: GTTGGTAAATCTCGTGCCAGC    
MiFish-E 12S amplicon R: CATAGTGGGGTATCTAATCCTAGTTTG    
  
#### Metadata sheet (optional) 

The metadata file has to follow the QIIME2 specifications (https://docs.qiime2.org/2021.2/tutorials/metadata/). Below is a preview of the sample sheet used for this test. Keep the column headers the same for future use. The first column needs to be "ID" and can only contain numbers, letters, or "-". This is different than the sample sheet. NAs should be empty cells rather than "NA". 

### Create samplesheet sheet for ampliseq 

This file indicates the sample ID and the path to R1 and R2 files. Below is a preview of the sample sheet used in this test. File created on RStudio Interactive on Discovery Cluster using (`create_metadatasheets.R`).  

- sampleID (required): Unique sample IDs, must start with a letter, and can only contain letters, numbers or underscores (no hyphons!).  
- forwardReads (required): Paths to (forward) reads zipped FastQ files  
- reverseReads (optional): Paths to reverse reads zipped FastQ files, required if the data is paired-end  
- run (optional): If the data was produced by multiple sequencing runs, any string  

*This is an R script, not slurm script. Open RStudio interactive on Discovery Cluster to run this script.*

Prior to running R script, use the `rawdata` file created for the fastqc slurm array from within the raw data folder to create a list of files. Below is an example from our Offshore Wind project but the specifics of the sampleID will be project dependent. This project had four sequencing runs with different file names. 

`01a-metadata.R`

```
## Creating sample sheet for offshore wind eDNA project 

### Step 1: In terminal 

# cd raw_data 
# ls * > ../metadata/sample_list.txt

### Resume steps below

library(dplyr)
library(stringr)
library(strex)
#library(filesstrings)

### Read in sample sheet 

sample_list <- read.delim2("/work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/rawdata", header=F) %>% 
  dplyr::rename(forwardReads = V1) %>%
  mutate(sampleID = str_after_nth(forwardReads, "data/", 1),
         sampleID = str_before_nth(sampleID, "_R", 1),
         sampleID = gsub("Degen", "", sampleID),
         sampleID = gsub("_L001", "", sampleID),
         sampleID = gsub("Bottom", "_B", sampleID),
         sampleID = gsub("Surface", "_S", sampleID),
         sampleID = gsub("NA", "_NA", sampleID),
         sampleID = gsub("2B", "2B_NA", sampleID),
         sampleID = gsub("2A", "2A_NA", sampleID),
         sampleID = gsub("1A", "1A_NA", sampleID),
         sampleID = gsub("1B", "1B_NA", sampleID),
         sampleID = gsub("Blank_1", "BK1", sampleID),
         sampleID = gsub("Blank_2", "BK2", sampleID),
         sampleID = ifelse(!grepl('July', sampleID), sub("_.*", "", sampleID), sampleID)
         )

# creating sample ID 
sample_list$sampleID <- gsub("-", "_", sample_list$sampleID)

# keeping only rows with R1
sample_list <- filter(sample_list, grepl("R1", forwardReads, ignore.case = TRUE))

# duplicating column 
sample_list$reverseReads <- sample_list$forwardReads

# replacing R1 with R2 in only one column 
sample_list$reverseReads <- gsub("R1", "R2", sample_list$reverseReads)

# rearranging columns 
sample_list <- sample_list[,c(2,1,3)]

sample_list %>% write.csv("/work/gmgi/Fisheries/eDNA/offshore_wind2023/metadata/samplesheet.csv", 
                          row.names=FALSE, quote = FALSE)
```

### Run nf-core/ampliseq (Cutadapt & DADA2)

Some parameters will change based on Riaz or MiFish primer set. 

`01b-ampliseq.sh`:

```
#!/bin/bash
#SBATCH --error=output_messages/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=ampliseq
#SBATCH --mem=70GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project 
## 2. Adjust SBATCH options above (time, mem, ntasks, etc.) as desired  
## 3. Fill in F and R primer information (no reverse compliment)
## 4. Adjust parameters as needed (below is Fisheries team default for 12S)

# LOAD MODULES
# module load singularity/3.10.3
# module load nextflow/23.10.1

# Activate conda environment
source ~/../../work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

# SET PATHS 
metadata="" 
output_dir=""

nextflow run nf-core/ampliseq -resume \
   -profile singularity \
   --input ${metadata}/samplesheet.csv \
   --FW_primer "" \
   --RV_primer "" \
   --outdir ${output_dir} \
   --trunclenf 100 \
   --trunclenr 100 \
   --trunc_qmin 25 \
   --max_len 200 \
   --max_ee 2 \
   --min_len_asv 100 \
   --max_len_asv 115 \
   --sample_inference pseudo \
   --skip_taxonomy \
   --ignore_failed_trimming
```

To run:   
- `sbatch 01b-ampliseq.sh` 

## Step 5: Blast ASV sequences (output from DADA2) against our 3 databases 

### Populating /work/gmgi/databases folder 

We use NCBI, Mitofish, and GMGI-12S databases. 

#### Download NBCI

**Option 1**: Download ncbi-blast+ to `/work/gmgi/packages` using `wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz` and then `tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz`. NCBI latest: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. Once downloaded, user does not need to repeat this. I was struggling with the remote flag here within a slurm script.

**Option 2**: Download 12S sequences from NCBI via CRABS: Creating Reference databases for Amplicon-Based Sequencing.

CRABS requires python 3.9 so I created a new conda environment with this version: `conda create --name Env_py3.9 python=3.9`. Fisheries eDNA uses python 3.12. I was struggling to download CRABS in either conda environment... come back to this.

https://github.com/gjeunen/reference_database_creator

```
conda activate Env_py3.9
conda install -c bioconda crabs

cd /work/gmgi/databases/12S/ncbi

crabs db_download --source ncbi --database nucleotide --query '12S[All Fields]' --output 12S_ncbi_[date].fasta --keep_original no --batchsize 5000

makeblastdb -in 12S_ncbi_[date].fasta -dbtype nucl -out 12S_ncbi_[date] -parse_seqids
```

**Option 3**: `blast` package is downloaded in the fisheries_eDNA conda environment with `conda install blast`. Install nt database `update_blastdb.pl --decompress nt` once inside the `work/gmgi/databases/ncbi/nt` folder. This is not ideal because it will take up more space and need to be updated.


#### Download Mitofish 

Check Mitofish webpage (https://mitofish.aori.u-tokyo.ac.jp/download/) for the most recent database version number. Compare to the `work/gmgi/databases/12S/reference_fasta/12S/Mitofish/` folder. If needed, update Mitofish database:

```
## download db 
wget https://mitofish.aori.u-tokyo.ac.jp/species/detail/download/?filename=download%2F/complete_partial_mitogenomes.zip  

## unzip 
unzip 'index.html?filename=download%2F%2Fcomplete_partial_mitogenomes.zip'

## clean headers 
awk '/^>/ {print $1} !/^>/ {print}' mito-all > Mitofish_v4.02.fasta

## remove excess files 
rm mito-all* 
rm index*

## make NCBI db 
## make sure fisheries_eDNA conda environment is activated or module load ncbi-blast+/2.13.0
makeblastdb -in Mitofish_v4.02.fasta -dbtype nucl -out Mitofish_v4.02 -parse_seqids
```

Alternate option: Download Mitofish db with CRABS. This program and will download and format the db accordingly.   

```
git clone https://github.com/gjeunen/reference_database_creator.git

## Download Mitofish 
crabs db_download --source mitofish --output /work/gmgi/databases/12S/Mitofish/mitofish.fasta --keep_original yes
### I couldn't get the function crabs to work 
```

#### Download GMGI 12S 




### Running taxonomic ID script 

`02-taxonomicID.sh`: 

```
#!/bin/bash
#SBATCH --error=output_messages/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=tax_ID
#SBATCH --mem=30GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project; change db path if not 12S

## LOAD MODULES 
## can use module on NU cluster or own ncbi-blast+
# module load ncbi-blast+/2.13.0
ncbi_program="/work/gmgi/packages/ncbi-blast-2.16.0+"

# Activate conda environment
source ~/../../work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

# SET PATHS 
ASV_fasta=""
out=""

gmgi="/work/gmgi/databases/12S/GMGI"
mito="/work/gmgi/databases/12S/Mitofish"
taxonkit="/work/gmgi/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
${ncbi_program}/bin/blastn -remote -db nt \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_NCBI.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore' \
   -verbose

## Mitofish database 

${ncbi_program}/bin/blastn -db ${mito}/*.fasta \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_Mito.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid  pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## GMGI database 
${ncbi_program}/bin/blastn -db ${gmgi}/*.fasta \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_GMGI.txt \
   -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

############################

#### TAXONOMIC CLASSIFICATION #### 
## creating list of staxids from all three files 
awk -F $'\t' '{ print $4}' ${out}/BLASTResults_NCBI.txt | sort -u > ${out}/NCBI_sp.txt

## annotating taxid with full taxonomic classification
cat ${out}/NCBI_sp.txt | ${taxonkit}/taxonkit reformat -I 1 -r "Unassigned" > ${out}/NCBI_taxassigned.txt
```

To run:  
- `sbatch 02-taxonomicID.sh` 