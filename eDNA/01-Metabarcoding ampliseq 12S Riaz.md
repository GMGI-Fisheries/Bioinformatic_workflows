# Metabarcoding workflow for 12S amplicon sequencing with Riaz primers

*page details in progress.* 

The 12S rRNA gene region of the mitogenome is ~950 bp. There are two popular primer sets to amplify two different regions of 12S: Riaz and MiFish. The following workflow includes script specific to the **Riaz** primer set.

![](https://th.bing.com/th/id/OIP.EbXPYETYLBPEymNEIVEGLQHaCc?rs=1&pid=ImgDetMain)

Citation: [Riaz et al. 2011](https://academic.oup.com/nar/article/39/21/e145/1105558)

Scripts to run: 

1. 00-fastqc.sh   
2. 00-multiqc.sh  
3. 01a-metadata.R
4. 01b-ampliseq.sh
5. 02-taxonomicID.sh  

## Step 1: Conda environment: Fisheries eDNA 

Background information on Conda: https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html. 

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
conda install blast
conda install singularity
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

# Update nextflow ampliseq workflow 
nextflow pull nf-core/ampliseq
``` 
 
## Step 2: Assess quality of raw data  

Background information on FASTQC: https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/05_qc_running_fastqc_interactively.html. 

`00-fastqc.sh`: 

```
#!/bin/bash
#SBATCH --error=output/fastqc_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/fastqc_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
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
- Start slurm array (e.g., with 138 files) = `sbatch --array=0-137 00-fastqc.sh`.

Notes:  
- This is going to output *many* error and output files. After job completes, use `cat *output.* > ../fastqc_output.txt` to create one file with all the output and `cat *error.* > ../fastqc_error.txt` to create one file with all of the error message outputs. 
- Within the `out_dir` output folder, use `ls *html | wc` to count the number of html output files (1st/2nd column values). This should be equal to the --array range used and the number of raw data files. If not, the script missed some input files so address this before moving on.  


## Step 3: Visualize quality of raw data  

Background information on MULTIQC: https://multiqc.info/docs/#:~:text=MultiQC%20is%20a%20reporting%20tool%20that%20parses%20results,experiments%20containing%20multiple%20samples%20and%20multiple%20analysis%20steps.

`00-multiqc.sh` 

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
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

[Nf-core](https://nf-co.re/): A community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/).  
Nextflow: scalable and reproducible scientific workflows using software containers, used to build wrapper programs like the one we use here.  

https://nf-co.re/ampliseq/2.11.0: nfcore/ampliseq is a bioinformatics analysis pipeline used for amplicon sequencing, supporting denoising of any amplicon and supports a variety of taxonomic databases for taxonomic assignment including 16S, ITS, CO1 and 18S. 

![](https://raw.githubusercontent.com/nf-core/ampliseq/2.11.0//docs/images/ampliseq_workflow.png)

We use ampliseq for the following programs:  
- [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200) is trimming primer sequences from sequencing reads. Primer sequences are non-biological sequences that often introduce point mutations that do not reflect sample sequences. This is especially true for degenerated PCR primer. If primer trimming would be omitted, artifactual amplicon sequence variants might be computed by the denoising tool or sequences might be lost due to become labelled as PCR chimera.  
- [DADA2](https://www.nature.com/articles/nmeth.3869) performs fast and accurate sample inference from amplicon data with single-nucleotide resolution. It infers exact amplicon sequence variants (ASVs) from amplicon data with fewer false positives than many other methods while maintaining high sensitivity.  

We skip the taxonomic assignment because we use 3-db approach described in the next section. 

### 12S primer sequences (required)

Below is what we used for 12S amplicon sequencing. Ampliseq will automatically calculate the reverse compliment and include this for us.

Riaz 12S amplicon F Original: ACTGGGATTAGATACCCC  
Riaz 12S amplicon F Degenerate: ACTGGGATTAGATACCCY     
Riaz 12S amplicon R: TAGAACAGGCTCCTCTAG     

### Metadata sheet (optional) 

The metadata file has to follow the QIIME2 specifications (https://docs.qiime2.org/2021.2/tutorials/metadata/). Below is a preview of the sample sheet used for this test. Keep the column headers the same for future use. The first column needs to be "ID" and can only contain numbers, letters, or "-". This is different than the sample sheet. NAs should be empty cells rather than "NA". 

### Create samplesheet sheet for ampliseq 

This file indicates the sample ID and the path to R1 and R2 files. Below is a preview of the sample sheet used in this test. File created on RStudio Interactive on Discovery Cluster using (`create_metadatasheets.R`).  

- sampleID (required): Unique sample IDs, must start with a letter, and can only contain letters, numbers or underscores (no hyphons!).  
- forwardReads (required): Paths to (forward) reads zipped FastQ files  
- reverseReads (optional): Paths to reverse reads zipped FastQ files, required if the data is paired-end  
- run (optional): If the data was produced by multiple sequencing runs, any string  

| sampleID | forwardReads              | reverseReads              | run |
|----------|---------------------------|---------------------------|-----|
| sample1  | ./data/S1_R1_001.fastq.gz | ./data/S1_R2_001.fastq.gz | A   |
| sample2  | ./data/S2_fw.fastq.gz     | ./data/S2_rv.fastq.gz     | A   |
| sample3  | ./S4x.fastq.gz            | ./S4y.fastq.gz            | B   |
| sample4  | ./a.fastq.gz              | ./b.fastq.gz              | B   |

*This is an R script, not slurm script. Open RStudio interactive on Discovery Cluster to run this script.*

Prior to running R script, use the `rawdata` file created for the fastqc slurm array from within the raw data folder to create a list of files. Below is an example from our Offshore Wind project but the specifics of the sampleID will be project dependent. This project had four sequencing runs with different file names. 

`01a-metadata.R`

```
## Load libraries 

library(dplyr)
library(stringr)
library(strex) 

### Read in sample sheet 

sample_list <- read.delim2("/work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/rawdata", header=F) %>% 
  dplyr::rename(forwardReads = V1) %>%
  mutate(sampleID = str_after_nth(forwardReads, "data/", 1),
         sampleID = str_before_nth(sampleID, "_R", 1))

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

Update ampliseq workflow if needed: `nextflow pull nf-core/ampliseq`. 

`01b-ampliseq.sh`:

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
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
## 3. Fill in F primer information based on primer type (no reverse compliment needed)
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
   --RV_primer "TAGAACAGGCTCCTCTAG" \
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

#### Files generated by ampliseq 

Pipeline summary reports:  
- `summary_report/`
- `summary_report.html`: pipeline summary report as standalone HTML file that can be viewed in your web browser.
- `*.svg*`: plots that were produced for (and are included in) the report.
- `versions.yml`: software versions used to produce this report.

Preprocessing:  
- FastQC: `fastqc/` and `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.  
- Cutadapt: `cutadapt/` and `cutadapt_summary.tsv`: summary of read numbers that pass cutadapt  
- MultiQC: `multiqc`, `multiqc_data/`, `multiqc_plots/` with `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser; 

ASV inferrence with DADA2:  
- `dada2/`, `dada2/args/`, `data2/log/` 
   - `ASV_seqs.fasta`: Fasta file with ASV sequences.
   - `ASV_table.tsv`: Counts for each ASV sequence.
   - `DADA2_stats.tsv`: Tracking read numbers through DADA2 processing steps, for each sample.
   - `DADA2_table.rds`: DADA2 ASV table as R object.
   - `DADA2_table.tsv`: DADA2 ASV table.  
- `dada2/QC/`
   - `*.err.convergence.txt`: Convergence values for DADA2's dada command, should reduce over several magnitudes and approaching 0.  
   - `*.err.pdf`: Estimated error rates for each possible transition. The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. The estimated error rates (black line) should be a good fit to the observed rates (points), and the error rates should drop with increased quality.  
   - `*_qual_stats.pdf`: Overall read quality profiles: heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position.  
   - `*_preprocessed_qual_stats.pdf`: Same as above, but after preprocessing.  

We add an ASV length filter that will output `asv_length_filter/` with:  
- `ASV_seqs.len.fasta`: Fasta file with filtered ASV sequences.  
- `ASV_table.len.tsv`: Counts for each filtered ASV sequence.  
- `ASV_len_orig.tsv`: ASV length distribution before filtering.  
- `ASV_len_filt.tsv`: ASV length distribution after filtering.  
- `stats.len.tsv`: Tracking read numbers through filtering, for each sample.  

## Step 5: Blast ASV sequences (output from DADA2) against our 3 databases 

### Populating /work/gmgi/databases folder 

We use NCBI, Mitofish, and GMGI-12S databases. 

#### Download and/or update NBCI blast nt database

NCBI is updated daily and therefore needs to be updated each time a project is analyzed. This is the not the most ideal method but we were struggling to get the `-remote` flag to work within slurm because I don't think NU slurm is connected to the internet? NU help desk was helping for awhile but we didn't get anywhere.

Within `/work/gmgi/databases/ncbi`, there is a `update_nt.sh` script with the following code. To run `sbatch update_nt.sh`. This won't take long as it will check for updates rather than re-downloading every time. 

```
#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=update_ncbi_nt
#SBATCH --mem=50G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate fisheries_eDNA

# Create output directory if it doesn't exist
cd /work/gmgi/databases/ncbi/nt

# Update BLAST nt database
update_blastdb.pl --decompress nt

# Print completion message
echo "BLAST nt database update completed"
```

View the `update_ncbi_nt.out` file to confirm the echo printed at the end.

*Emma is still troubleshooting the -remote flag to also avoid storing the nt db within our /work/gmgi folder.* 

#### Download and/or update Mitofish database  

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
makeblastdb -in Mitofish_v4.02.fasta -dbtype nucl -out Mitofish_v4.02.fasta -parse_seqids
```

#### Download GMGI 12S 

This is our in-house GMGI database that will include version numbers. Check `/work/gmgi/databases/12S/GMGI/` for current uploaded version number and check our Box folder for the most recent version number. 

On OOD portal, click the Interactive Apps dropdown. Select Home Directory under the HTML Viewer section. Navigate to the `/work/gmgi/databases/12S/GMGI/` folder. In the top right hand corner of the portal, select Upload and add the most recent .fasta file from our Box folder. 

To create a blast db from this reference fasta file (if updated): 

```
cd /work/gmgi/databases/12S/GMGI/ 

## make NCBI db 
## make sure fisheries_eDNA conda environment is activated 
### CHANGE THE VERSION NUMBER BELOW TO LATEST
makeblastdb -in GMGI_Vert_Ref_2024v1.fasta -dbtype nucl -out GMGI_Vert_Ref_2024v1.fasta
```

### Running taxonomic ID script 

`02-taxonomicID.sh`: 

```
#!/bin/bash
#SBATCH --error=output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=tax_ID
#SBATCH --mem=30GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

### USER TO-DO ### 
## 1. Set paths for project; change db path if not 12S

# Activate conda environment
source /work/gmgi/miniconda3/bin/activate
conda activate fisheries_eDNA

# SET PATHS 
ASV_fasta=""
out=""

gmgi="/work/gmgi/databases/12S/GMGI"
mito="/work/gmgi/databases/12S/Mitofish"
ncbi="/work/gmgi/databases/ncbi/nt"
taxonkit="/work/gmgi/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
blastn -db ${ncbi}/*.fasta \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_NCBI.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore' \
   -verbose

## Mitofish database 
blastn -db ${mito}/*.fasta \
   -query ${ASV_fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_Mito.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid  pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## GMGI database 
blastn -db ${gmgi}/*.fasta \
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