# Base scripts for metabarcoding 

Scripts are located on NU's Discovery Cluster in `/work/gmgi/Fisheries/base_scripts/metabarcoding`. 

*page details in progress.* 

Scripts to run: 

1. 00-fastqc.sh   
2. 00-multiqc.sh  
3. 01a-metadata.R
4. 01b-ampliseq.sh
5. 02-taxonomicID.sh  

## Copy script folder to project folder 

This will copy all scripts in the `base_scripts/metabarcoding` folder into your specific project folder. Replace `project_path` with the path to specific project folder. 

`cp /work/gmgi/Fisheries/base_scripts/metabarcoding/* project_path`

## Contents 

### Step 1: Assess quality of raw data  

`00-fastqc.sh`: 

```
#!/bin/bash
#SBATCH --error=output_messages/fastqc_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/fastqc_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=fastqc
#SBATCH --mem=3GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for your project

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


### Step 2: Visualize quality of raw data  

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

## SET PATHS 
## fastqc_output = output from 00-fastqc.sh; fastqc program
fastqc_output=""
multiqc_dir=""

## RUN MULTIQC 
multiqc --interactive ${fastqc_output} -o ${multiqc_dir} --filename multiqc_raw.html
```

To run:  
- Activate conda environment `conda activate haddock_methylation` (to be changed; use haddock_methylation for now)  
- `sbatch 00-multiqc.sh` 

Notes:  
- Depending on the number of files per project, multiqc can be quick to run without a slurm script. To do this, run each line separately in the command line after activating the conda environment.  

### Step 3: Create metadata sheet for ampliseq 

*This is an R script, not slurm script. Open RStudio interactive on Discovery Cluster to run this script.*

Prior to running R script, use `ls * > ../metadata/sample_list.txt` from within the raw data folder to create a list of files.

`01a-metadata.R`

```
## Creating sample sheet for ampliseq

### Step 1: In terminal 

# cd raw_data 
# ls * > ../metadata/sample_list.txt

### Resume steps below

library(dplyr)
library(stringr)
library(strex)
#library(filesstrings)

### USER TO-DO ### 
## 1. Complete Step 1 above to create file name list 
## 2. Set paths for your project 
## 3. Add other sampleID manipulation as needed 
## 4. Edit metadata path in export.csv command

### Set paths (USER EDITS)
sample_list_path=""
raw_data_path=""

### Read in sample sheet 
sample_list <- read.delim2(sample_list_path, header=F) %>% 
  dplyr::rename(forwardReads = V1) %>% 
  filter(!forwardReads == "sample_list.txt") %>%
  filter(!forwardReads == "rawdata")

### creating sample ID 
sample_list$sampleID <- gsub("-", "_", sample_list$forwardReads)

### Adding other sampleID manipulation (USER EDITS)
### examples below 
# sample_list$sampleID <- str_before_nth(sample_list$forwardReads, "-12S", 1)
# sample_list <- sample_list %>%
#  mutate(sampleID = ifelse(grepl('H2O', forwardReads), "H2O-negative", sampleID))

### keeping only rows with R1
sample_list <- filter(sample_list, grepl("R1", forwardReads, ignore.case = TRUE))

### duplicating column 
sample_list$reverseReads <- sample_list$forwardReads

### replacing R1 with R2 in only one column 
sample_list$reverseReads <- gsub("R1", "R2", sample_list$reverseReads)

# rearranging columns 
sample_list <- sample_list[,c(2,1,3)]

# adding file path
sample_list$forwardReads <- paste(raw_data_path, sample_list$forwardReads, sep = "")
sample_list$reverseReads <- paste(raw_data_path, sample_list$reverseReads, sep = "")

## exporting as csv (USER EDITS PATH PRIOR TO METADATA)
sample_list %>% write.csv("/metadata/samplesheet.csv", row.names=FALSE, quote = FALSE)
```

### Step 4: Run nf-core/ampliseq (Cutadapt & DADA2)

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
module load singularity/3.10.3
module load nextflow/23.10.1

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


### Step 5: Blast ASV sequences (output from DADA2) against our 3 databases 

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
module load ncbi-blast+/2.13.0

# SET PATHS 
ASV_fasta=""
out=""
db="/work/gmgi/databases/12S/reference_fasta"
taxonkit="/work/gmgi/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
blastn -remote -db nt \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_NCBI.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## Mitofish database 

blastn -db ${db}/Mitofish.fasta \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_Mito.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid  pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## GMGI database 
blastn -db ${db}/GMGIVertRef.fasta \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_GMGI.txt \
   -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

############################

#### TAXONOMIC CLASSIFICATION #### 
## creating list of staxids from all three files 
awk -F $'\t' '{ print $4}' ${out}/BLASTResults_NCBI.txt | sort -u > ${out}/NCBI_sp.txt

## annotating taxid with full taxonomic classification
cat ${out}/NCBI_sp.txt | ${taxonkit}/taxonkit reformat -I 1 -r "Unassigned" > NCBI_taxassigned.txt
```

To run:  
- `sbatch 02-taxonomicID.sh` 

Notes:  
- Emma currently troubleshooting running ncbi -remote within a slurm script. Can run these lines of code just on command line and will work fine.
