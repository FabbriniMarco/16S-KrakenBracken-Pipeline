<div align=center>
  <img width="196" alt="Circle logo" src="https://github.com/user-attachments/assets/7633e817-fb7f-4b11-97d6-530d808275c6" height="186">
</div>


# 16S-KrakenBracken-Pipeline

This repository hosts a [Nextflow](https://www.nextflow.io/) containerized pipeline tailored for 16S rRNA amplicon sequencing analyses, leveraging [Kraken2](https://github.com/DerrickWood/kraken2)/[Bracken](https://github.com/jenniferlu717/Bracken) to generate abundance tables (OTU tables) and compute alpha and beta diversity results.

The pipeline is designed to process FASTQ 16S amplicon files, and it can handle mixed paired-end and single-end data. Upon input, sequences undergo quality filtering using [fastp](https://github.com/OpenGene/fastp), ensuring high-quality data for downstream analysis. Additionally, the pipeline offers the capability to filter sequences for a reference genome using [KneadData](https://github.com/biobakery/kneaddata), providing flexibility for various research needs.

Processed sequences are then mapped onto the [SILVA 138_1 SSURef database](https://www.arb-silva.de/no_cache/download/archive/release_138_1/) using Kraken2, with taxonomic assignments obtained through Bracken. The resulting taxonomic outputs are further processed to generate counts and relative percentages for OTU tables. Additionally, the pipeline calculates alpha diversity values using several metrics and produces beta diversity distance matrices, enriching the depth of analysis. The outputs are fully compatible and ready to use with the [microbAIDeR R-package](https://github.com/FabbriniMarco/microbAIDeR) that you can find among my repos.
A schematic overview of the pipeline is shown in the chart below, created in Lucidchart ([Lucid.co](https://lucid.co/)):


<div align=center>
  <img alt="Pipeline schematics" src="https://github.com/user-attachments/assets/4fc4766f-8102-4268-9812-51c58048c89c" >
</div>

<br>

## Table of contents
1. Using Nextflow (suggested)
	1. [Installing Nextflow](#installing-nextflow)
 	2. [Obtaining the pipeline](#obtaining-the-pipeline)
  	3. [Launch nextflow pipeline](#launch-nextflow-pipeline)
   	4. [Further customization options](#further-customization-options)
   	5. [Generating a bowtie2-compatible reference genome](#generating-a-bowtie2-compatible-reference-genome)
   	6. [Nextflow pipeline outputs](#nextflow-pipeline-outputs)
2. [Use the legacy script from the container](#use-the-legacy-script-from-the-container)
	1. [Dependencies](#dependencies)
	2. [Installation](#installation)
		- [Pre-built container (recommended)](#pre-built-container-recommended)
		- [Build your own container](#build-your-own-container)
	3. [Usage instructions](#usage-instructions)
		- [Running inside the container](#running-inside-the-container)
		- [Pass a script to the container](#pass-a-script-to-the-container)
	4. [Legacy pipeline outputs](#legacy-pipeline-outputs)
	5. [Generating a BMTagger-compatible reference genome](#generating-a-bmtagger-compatible-reference-genome)
	6. [Parameters of the legacy script](#parameters-of-the-script)
3. [Citation](#citation)

<br>

***

## Using Nextflow

### Installing Nextflow
Starting from version 0.2 of this pipeline, Nextflow can be used to quickly and smoothly run the pipeline on your machine. Of course, Nextflow needs to be installed in the machine. You can consult Nextflow Documentation [here](https://www.nextflow.io/docs/latest/install.html) for installation instruction. <br>
Briefly, if you don't have Java 11+ installed on your system you can install it via SDKMAN as such:
```bash
curl -s https://get.sdkman.io | bash
sdk install java 17.0.10-tem
# Confirm java installation
java -version
```
And then install nextflow:
```bash
# This command will create a nextflow execution in the current folder
curl -s https://get.nextflow.io | bash
chmod +x nextflow
# Move nextflow executable to a folder in $PATH, for example:
mv nextflow ~/.local/bin/
# or for a system-wise installation (requires sudo privileges)
mv nextflow /usr/bin/
```

<br>

### Obtaining the pipeline
In order to use the pipeline you need to have two files: the ```nextflow.config``` and ```kraken16s.nf```. You can obtain them in several ways:
- download them from the [Release](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/releases) section of this repo clicking [here](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/releases/download/v0.2/16s-krakenbracken-pipeline-v2.tar.gz). Extract the archive (```tar -zxvf 16s-krakenbracken-pipeline-v2.tar.gz```)
- clone this repository with ```git clone https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline```
- manually download the nextflow [main file](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/blob/main/kraken16S.nf) and [config file](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/blob/main/nextflow.config) from this repo

<br>

Once you have obtained the files you can add the folder to $PATH to be able to launch the pipeline from anywhere on your machine:
```bash
echo 'export PATH="/path/to/16s-krakenbracken-pipeline-v2/:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

To test the installation you can run
```
nextflow run kraken16s.nf --help
```
You should see the help message:
<div align=center>
  <img alt="Help section" src="https://github.com/user-attachments/assets/96d36e50-836f-4811-b398-ce6aa3e96b48" >
</div>

<br>

### Launch nextflow pipeline
After consulting the help section, you can start the pipeline with default steps (_i.e._ performing both ```fastp``` read processing and host genome filtering with ```KneadData```) providing an input folder and a reference genome path as such:
``` bash
nextflow run kraken16s.nf --fastq_folder raw_seqs --genome_path HOST_BOWTIE2_DB_PATH
```
Remember to specify in the ```--genome_path``` argument the path of a *bowtie2-formatted* host genome for the KneadData filtering step (See below for the guide on [How to generate a bowtie2-compatible reference genome](#generating-a-bowtie2-compatible-reference-genome).
<br>
<br>
If you don't want to perform host filtering, the process can be skipped by setting the use_kneaddata argument to false. Similarly you can skip fastp read processing, if you want:
```bash
# Run the pipeline without host genome filtering
nextflow run kraken16s.nf --fastq_folder raw_seqs --use_kneaddata false
# Run the pipeline without host genome filtering and without clipping reads with fastp
nextflow run kraken16s.nf --fastq_folder raw_seqs --use_kneaddata false --use_fastp false
```

If for any reason the pipeline gets stopped, you may want to resume it by adding the ```-resume``` flag in the nextflow run command for using previously cached contents.
<br>
<br>


### Further customization options
In the commands above, the pipeline uses all available cores on the machine, spawning all the processes possible. 
In order to control the resource usage, you can act on the nextflow call itself:
```bash
# Run the script using a maximum of 20 threads and 32GB of system RAM, spawning a maximum of 20 parallel instances
nextflow run kraken16s.nf --fastq_folder raw_seqs --genome_path "/mnt/databases/bowtie2_hg38" -process.cpus 20 -process.memory '32 GB' -process.maxForks 20
```

In addition, you can specify a custom kraken2-compatible database for 16S analyses, which has to be formatted for Bracken using a ```-l ${READ_LEN}``` parameter of 250 or 500. By default the pipeline uses Kraken2 special database [Silva 138.1 SSuRef NR99](https://www.arb-silva.de/). You can consult Kraken2 section concerning 16S databases [here](https://github.com/DerrickWood/kraken2/wiki/Manual#16s-databases).

<br>
<br>

### Generating a bowtie2-compatible reference genome
In order to filter the host genome using the nextflow pipeline you need to have a local copy of the host genome, indexed for bowtie2.
```bash
# Download a reference genome in FASTA format, for example the human genome
cd /mnt/databases # Custom directory of your choice
mkdir bowtie2_hg38 && cd bowtie2_hg38
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz

# Build the database using bowtie2-build
bowtie2-build -f --threads 40 GRCh38_latest_genomic.fna bowtie2_hg38

# Optionally remove the fasta file
rm -fr GRCh38_latest_genomic.fna

# You should see several bowtie2_hg38.*.bt2 files
pwd
# Use this path as 'HOST_BOWTIE2_DB_PATH' in the nextflow run command
```
You can have a look at the [bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) for further details.

<br>
<br>

### Nextflow pipeline outputs

The output of the Nextflow pipeline looks like this:

``` diff
AnalysisKraken_${currentDate}/
├── fastp   # Contains samples folder with gzipped reads after fastp processing (if use_fastp : 'true')
├── kraken2  # Contains samples folders with Kraken2 mapping .kreport files
├── bracken  # Contains samples folders with Bracken outputs at phylum, family and genera levels (L2, L5, L6)
├── RESULTS/
    ├── alpha.xlsx  # Contains several alpha diversity metrics value for each sample
    ├── beta_diversity/  # Each subfolder contains a pairwise distance matrix. Ready-to-use for post-processing with microbAIDeR compute_beta_diversity() R-function 
    │   ├── braycurtis/
    │   │   └── distance-matrix.tsv
    │   ├── jaccard/
    │   │   └── distance-matrix.tsv
    │   ├── unweighted_unifrac/
    │   │   └── distance-matrix.tsv
    │   └── weighted_unifrac/
    │       └── distance-matrix.tsv
    ├── final_otu_table_matrix.csv  # Otu table with counts reported across levels L2, L5 and L6: (L2, phylum; L5, family; L6, genera)
    ├── final_otu_table_matrix_relabb.csv  # Otu table with relative counts reported across levels L2, L5 and L6. The relative abundance sums up to 100 within each level individually
    ├── otu_table_final.tsv  # Otu table as "final_otu_table_matrix.csv" but in tidy format
    └── taxa_summary/  # Otu table in tabular format divided by taxonomic ranks
        ├── otu_table_L2.tsv  
        ├── otu_table_L5.tsv
        └── otu_table_L6.tsv

```

<br>


After the nextflow run is completed, you can get rid of temporary files - if you don't need them - by running:
```bash
nextflow clean -f -q

# You can also delete the fastp, kraken2 and bracken folders
rm -fr AnalysisKraken_*/fastp
rm -fr AnalysisKraken_*/kraken2
rm -fr AnalysisKraken_*/bracken
```

<br>
<br>

***

## Use the legacy script from the container
Instead of using the Nextflow pipeline, you can also download the singularity container used by the pipeline and perform the analysis using the ```kraken16S.sh``` script instead. The logic of the workflow is the same; however, the legacy shell script performs host filtering using ```BMTagger``` with KneadData instead of bowtie2. In addition, the data output contains more files, which may not be necessary.
<br>

### Dependencies

Before using the container, ensure that you have [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) installed on your machine.

<br>

### Installation

#### Pre-built container (recommended)

The prebuilt ready-to-use container is [available on Sylabs cloud platform](https://cloud.sylabs.io/library/fabbrinimarco/16s-krakenbracken-pipeline/k16s) and can be easily downloaded: 

```bash
singularity pull library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:v2
# For an added layer of security, it can be also pulled by unique sha256 sum
singularity pull --arch amd64 library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:sha256.af116e7b1e7392b5742035d7226ca5ffbd5b40eefdb8f0ae3e3d66b142edfb23
# Check the container
singularity run k16s_v2.sif kraken16S 
```

If you want to download the original version of the container, you can still pull the ```k16s_v1.sif``` from SCS:
```bash
singularity pull --arch amd64 library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:v1
# or with shasum checks
singularity pull --arch amd64 library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:sha256.7dc4939abc3f620d47e06edeba9adea39d8975f2504ee88b74642a209910d917
```


<br>

#### Build your own container
Clone this repository or download the build_container.def file and use it to build the Singularity container:
```bash
singularity build k16s_v2.sif build_container.def
```
Note that during the Kraken2 Silva database building process, 24 threads will be utilized.
Due to the inconsistent access to NCBI FTP services, the container build might fail downloading NCBI BLAST+.
If you encounter this issue and still want to build your own container you can build it without BLAST+ as sandbox using the kraken16S_container_noBLAST+.def from this repository, install BLAST+ manually and then recreate the SIF image. Note that this approach will build the k16s_v1.sif container which is perfectly fine for using the legacy shell script. The difference between v1 and v2 containers is the presence of bowtie2 for the Nextflow pipeline. If you're in this section of the README, you are likely not interested in this feature.

```bash
singularity build --fix-perms --sandbox k16s_v1 kraken16S_container_noBLAST+.def
singularity exec --writable k16s_v1 bash
```

```bash
# Inside the container shell
cd /opt/programs

curl https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz.md5 -o known_checksum.md5
curl -o ncbi-blast-2.15.0+-x64-linux.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz

md5sum ncbi-blast-2.15.0+-x64-linux.tar.gz > curled_checksum.md5
diff known_checksum.md5 curled_checksum.md5
# If nothing prints, checksums correspond and you can proceed

rm -fr *.md5
tar -zxvf ncbi-blast-2.15.0+-x64-linux.tar.gz
# If the above command return EOF errors the file is deleted and should be downloaded once again

rm -fr ncbi-blast-2.15.0+-x64-linux.tar.gz
mv ncbi-blast-2.15.0+ ncbi-blast

exit
```

```bash
singularity build k16s_v1.sif k16s_v1/
```

<br>

### Usage instructions
#### Running inside the container
Executing commands within the container is straightforward. Assuming you have a folder named 'rawseqs' containing your raw 16S FASTQ sequences (in .fastq or .fastq.gz format) in the current path, run:

``` bash
singularity exec --bind $(pwd):/workdir k16s_v2.sif bash
cd /workdir
kraken16S -i rawseqs/
# Once the script is completed remember to close the container shell
exit
```

You can also perform reference genome filtering, such as with the human genome, by binding a local BMTagger-indexed genome to the container folder '/mnt/databases/bmtagger', as follows:

``` bash
# In this example I have a BMTagger-indexed database in the local path '/mnt/luks/databases/hg38'
singularity exec --bind /mnt/luks/databases/hg38:/mnt/databases/bmtagger --bind $(pwd):/workdir k16s_v2.sif bash
cd /workdir
kraken16S -i rawseqs/ -f TRUE
```
<br>

#### Pass a script to the container
You can also pass the container a shell script to be executed to avoid opening a shell inside the container:

```bash
echo -e "#!/bin/bash
cd /workdir
kraken16S -i rawseqs/" >> kraken16Scall.sh

chmod +x kraken16Scall.sh

singularity exec --bind $(pwd):/workdir k16s_v2.sif kraken16Scall.sh
```

Or with genome filtering, as the example above:

```bash
echo -e "#!/bin/bash
cd /workdir
kraken16S -i rawseqs/ -f TRUE" >> kraken16Scall.sh

chmod +x kraken16Scall.sh

singularity exec --bind /mnt/luks/databases/hg38:/mnt/databases/bmtagger --bind $(pwd):/workdir k16s_v2.sif kraken16Scall.sh
```

<br>

### Legacy pipeline outputs

The output of the pipeline looks like this:

``` diff
AnalysisKraken16_2024-05-13_silvaNR99/
├── paired_sample_list.txt
├── pipeline_log.txt  #Logs for debugging purposes
├── reads_counts_FILTERED.txt  # The number of reads in the forward FASTQ after FASTP filtering (if performed)
├── reads_counts.txt  # The number of reads in the RAW forward FASTQ
├── RESULTS
│   ├── alpha.xlsx  # Contains several alpha diversity metrics value for each sample
│   ├── beta_diversity  # Each subfolder contains a pairwise distance matrix
│   │   ├── braycurtis
│   │   │   └── distance-matrix.tsv
│   │   ├── jaccard
│   │   │   └── distance-matrix.tsv
│   │   ├── unweighted_unifrac
│   │   │   └── distance-matrix.tsv
│   │   └── weighted_unifrac
│   │       └── distance-matrix.tsv
│   ├── final_otu_table_matrix.csv  # Otu table with counts reported across levels L2, L5 and L6: (L2, phylum; L5, family; L6, genera)
│   ├── final_otu_table_matrix_relabb.csv  # Otu table with relative counts reported across levels L2, L5 and L6. The relative abundance sums up to 100 within each level individually
│   ├── logs.log  # Summary of the parameters passed to the call function
│   ├── otu_table_final.tsv  # Otu table in tidy format, as "final_otu_table_matrix.csv"
│   └── taxa_summary  # Otu table in tabular format divided by taxonomic ranks
│       ├── otu_table_L2.tsv  
│       ├── otu_table_L5.tsv
│       └── otu_table_L6.tsv
└── sample_list.txt
```

<br>

### Generating a BMTagger-compatible reference genome
To build a reference genome to be used for filtering the RAW reads with the legacy shell script, the database needs to be formatted properly.

```bash
mkdir hg38 && cd hg38
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
pigz -p 10 -d GRCh38_latest_genomic.fna.gz
mv GRCh38_latest_genomic.fna hg38.fa
bmtool -d hg38.fa -o hg38.bitmask -A 0 -w 18
srprism mkindex -i hg38.fa -o hg38.srprism -M 7168
makeblastdb -in hg38.fa -dbtype nucl
```
The database local path needs to be bound to the container in the /mnt/database/bmtagger path

<br>

### Parameters of the legacy script
You can consult the function details of the pipeline by simply invoking the main function name, with no further parameters required:

```bash
singularity exec k16s_v2.sif kraken16S
```


![Screenshot 2024-05-13 175759](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/478a84ea-c17f-4b21-a6b4-cab3a0d8f2dc)


For example, you can control the number of parallel threads used for the read mapping and read compress/decompress operations. You can decide whether or not to perform fastp read quality check and filtering, as well as eventual host genome filtering. If you deal with sequences with a different length than 250bp for paired-end data and 500bp for single-end data you can also edit the read length Bracken settings choosing either 250 or 500bp options.
You can also decide to perform only taxonomic assignment and skip the alpha and beta diversity analyses, or keep the temp files (i.e., the kraken2 and bracken outputs, as well as filtered reads).

<br>
<br>

***

## Citation

If you use this pipeline, please cite it as such:

```diff
@Manual{,
  title = {16S-KrakenBracken-Pipeline: Containerized pipeline for 16S rRNA amplicon sequencing analyses using Kraken2/Bracken},
  author = {Marco Fabbrini},
  year = {2024},
  note = {Singularity image v2},
  url = {https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline.git},
}
```
Example: 
> Fabbrini, M. (2024) '16S-KrakenBracken-Pipeline - Containerized pipeline for 16S rRNA amplicon sequencing analyses using Kraken2/Bracken'. Available on GitHub: https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline














