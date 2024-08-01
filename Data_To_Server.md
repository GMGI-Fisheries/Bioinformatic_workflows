# Downloading sequencing data to servers

## Downloading data from Illumina BaseSpace to HPC

https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview#InstallBaseSpaceSequenceHubCLI

### NU Discovery Cluster

1. Create folder called bin: `mkdir $HOME/bin`  
2. Download BaseSpace CLI: `wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O $HOME/bin/bs`  
3. Change the file permissions to make the downloaded binary executable: `chmod u+x $HOME/bin/bs`  
4. Authenticate your account: `bs auth`  
5. Navigate to the webpage provided and authenticate use of BaseSpace CLI.  
6. Find the Project ID of desired download: Within Sequence Hub, navigate to Projects. Select the desired project. The Project ID is in the webpage handle (e.g., https://basespace.illumina.com/projects/411088072/about). 
7. Navigate to `$HOME/bin` and download dataset: `bs download project -i <project_id> --extension=fastq.gz -o /local/output/path`.   
8. Navigate to the output path and move all files out of subdirectories: `mv */* .` 