# 16S-KrakenBracken-Pipeline

This repository hosts a containerized pipeline tailored for 16S rRNA amplicon sequencing analyses, leveraging [Kraken2](https://github.com/DerrickWood/kraken2)/[Bracken](https://github.com/jenniferlu717/Bracken) to generate abundance tables (OTU tables) and compute alpha and beta diversity results.

The pipeline is designed to process FASTQ 16S amplicon files, preferably in paired-end format. While it can handle mixed paired-end and single-end data, it does not support solely single-end data. Upon input, sequences undergo quality filtering using [fastp](https://github.com/OpenGene/fastp), ensuring high-quality data for downstream analysis. Additionally, the pipeline offers the capability to filter sequences for a reference genome using [KneadData](https://github.com/biobakery/kneaddata), providing flexibility for various research needs.

Processed sequences are then mapped onto the [SILVA 138_1 SSURef database](https://www.arb-silva.de/no_cache/download/archive/release_138_1/) using Kraken2, with taxonomic assignments obtained through Bracken. The resulting taxonomic outputs are further processed to generate counts and relative percentages for OTU tables. Additionally, the pipeline calculates alpha diversity values using several metrics and produces beta diversity distance matrices, enriching the depth of analysis.
A schematic overview of the pipeline is shown in the chart below, created in Lucidchart ([Lucid.co](https://lucid.co/)):

![image](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/ffa912cf-f19d-4320-8203-b7316bdc6075)

## Table of contents
1. [Dependencies](README.md#dependencies)
2. [Installation](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#installation)
	- [Pre-built container (recommended)](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#pre-built-container-recommended)
	- [Build your own container](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#build-your-own-container)
3. [Usage instructions](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#usage-instruction)
	- [Running inside the container](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#running-inside-the-container)
4. [Pipeline outputs](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#pipeline-outputs)
5. [Generating a BMTagger-compatible reference genome](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#generating-a-bmtagger-compatible-reference-genome)
6. [Parameters of the script](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline#parameters-of-the-script)

<br>

## Dependencies

Before using the container, ensure to have [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) installed on your machine.

<br>

## Installation

### Pre-built container (recommended)

The prebuilt ready-to-use container is [available on Sylabs cloud platform](https://cloud.sylabs.io/library/fabbrinimarco/16s-krakenbracken-pipeline/k16s) and can be easily downloaded: 

```bash
singularity pull library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:v1
# For an added layer of security, it can be also pulled by unique sha256 sum
singularity pull --arch amd64 library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:sha256.a5e72d21fe14fc4396ccd17117003f27e231f1c5af4559b1d53e8ba3415bed4c
# Check the container
singularity run K16S.sif kraken16S 
```

<br>

### Build your own container
Clone this repository or download the build_container.def file and use it to build the Singularity container:
```bash
singularity build K16S.sif build_container.def
```
Note that during the Kraken2 Silva database building process, 24 threads will be utilized.

<br>

## Usage instruction
### Running inside the container
Executing commands within the container is straightforward. Assuming you have a folder named 'rawseqs' containing your raw 16S FASTQ sequences (in .fastq or .fastq.gz format) in the current path, run:

``` bash
singularity exec --bind $(pwd):/workdir K16S.sif bash
cd /workdir
kraken16S -i rawseqs/
```

You can also perform reference genome filtering, such as with the human genome, by binding a local BMTagger-indexed genome to the container folder '/mnt/databases/bmtagger', as follows:

``` bash
# In this example I have a BMTagger-indexed database in the local path '/mnt/luks/databases/hg38'
singularity exec --bind /mnt/luks/databases/hg38:/mnt/databases/bmtagger --bind $(pwd):/workdir K16S.sif bash
cd /workdir
kraken16S -i test_data -f TRUE
```

<br>

## Pipeline outputs

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
│   ├── final_otu_table_matrix_relabb.csv  # Otu table with relative counts reported across levels L2, L5 and L6. The relativa abundance sums up to 100 within each level individually
│   ├── logs.log  # Summary of the parameters passed to the call function
│   ├── otu_table_final.tsv  # Otu table in tidy format, as "final_otu_table_matrix.csv"
│   └── taxa_summary  # Otu table in tabular format divided by taxonomic ranks
│       ├── otu_table_L2.tsv  
│       ├── otu_table_L5.tsv
│       └── otu_table_L6.tsv
└── sample_list.txt
```

<br>

## Generating a BMTagger-compatible reference genome
To build a reference genome to be used for filtering the RAW reads, the database needs to be formatted properly.

```bash
mkdir hg38 && cd hg38
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
pigz -p 10 -d GRCh38_latest_genomic.fna.gz
mv GRCh38_latest_genomic.fna hg38.fa
bmtool -d hg38.fa -o hg38.bitmask -A 0 -w 18
srprism mkindex -i hg38.fa -o hg38.srprism -M 7168
makeblastdb -in hg38.fa -dbtype nucl
```
The database local path needs to be binded to the container in the /mnt/database/bmtagger path

<br>

## Parameters of the script
You can consult the function details of the pipeline invoking the main function name:

```bash
singularity exec K16S.sif kraken16S
```
![image](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/523df061-48a8-4a38-9907-ffdfa5b682d9)














