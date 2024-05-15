# 16S-KrakenBracken-Pipeline

This repository hosts a containerized pipeline tailored for 16S rRNA amplicon sequencing analyses, leveraging [Kraken2](https://github.com/DerrickWood/kraken2)/[Bracken](https://github.com/jenniferlu717/Bracken) to generate abundance tables (OTU tables) and compute alpha and beta diversity results.

The pipeline is designed to process FASTQ 16S amplicon files, preferably in paired-end format. While it can handle mixed paired-end and single-end data, it does not support solely single-end data. Upon input, sequences undergo quality filtering using [fastp](https://github.com/OpenGene/fastp), ensuring high-quality data for downstream analysis. Additionally, the pipeline offers the capability to filter sequences for a reference genome using [KneadData](https://github.com/biobakery/kneaddata), providing flexibility for various research needs.

Processed sequences are then mapped onto the [SILVA 138_1 SSURef database](https://www.arb-silva.de/no_cache/download/archive/release_138_1/) using Kraken2, with taxonomic assignments obtained through Bracken. The resulting taxonomic outputs are further processed to generate counts and relative percentages for OTU tables. Additionally, the pipeline calculates alpha diversity values using several metrics and produces beta diversity distance matrices, enriching the depth of analysis.
A schematic overview of the pipeline is shown in the chart below, created in Lucidchart ([Lucid.co](https://lucid.co/)):

![image](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/ffa912cf-f19d-4320-8203-b7316bdc6075)

<br>

## Dependencies

Before using the container, ensure to have [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) installed on your machine.
<br><br>

## Installation

### Pre-built container (recommended)

Download the latest release from the [Release](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/releases) section. Extract the downloaded archive:

```bash
# Example
wget https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/releases/download/v1/K16S.sif.gz
pigz -d K16S.sif.gz
# Or
gunzip K16S.sif.gz
```
You will find the container inside the tarball. You can also download the container from the file list in the main branch of this repo.
<br>
### Build your own container
Clone this repository or download the build_container.def file and use it to build the Singularity container:
```bash
singularity build K16S.sif build_container.def
```
Note that during the Kraken2 Silva database building process, 24 threads will be utilized.
<br><br>
## Usage instruction
### Running inside the container
Executing commands within the container is straightforward. Assuming you have a folder named 'rawseqs' containing your raw 16S FASTQ sequences (in .fastq or .fastq.gz format) in the current path, run:

``` bash
singularity exec --bind $(pwd):/workdir K16S.sif bash
cd /workdir
kraken16S -i rawseqs/
```

You can also perform reference genome filtering, such as with the human genome, by binding a BMTagger-indexed genome located in the folder '/mnt/databases/bmtagger', as follows:

``` bash
# In this example I have a BMTagger-indexed database in the local path '/mnt/luks/databases/hg38'
singularity exec --bind /mnt/luks/databases/hg38:/mnt/databases/bmtagger --bind $(pwd):/workdir K16S.sif bash
cd /workdir
kraken16S -i test_data -f TRUE
```
<br><br>
## Pipeline outputs

The output of the pipeline looks like this:

``` diff
AnalysisKraken16_2024-05-13_silvaNR99/
├── paired_sample_list.txt
├── pipeline_log.txt
├── reads_counts_FILTERED.txt
├── reads_counts.txt
├── RESULTS
│   ├── alpha.xlsx
│   ├── beta_diversity
│   │   ├── braycurtis
│   │   │   └── distance-matrix.tsv
│   │   ├── jaccard
│   │   │   └── distance-matrix.tsv
│   │   ├── unweighted_unifrac
│   │   │   └── distance-matrix.tsv
│   │   └── weighted_unifrac
│   │       └── distance-matrix.tsv
│   ├── final_otu_table_matrix.csv
│   ├── final_otu_table_matrix_relabb.csv
│   ├── logs.log
│   ├── otu_table_final.tsv
│   └── taxa_summary
│       ├── otu_table_L2.tsv
│       ├── otu_table_L5.tsv
│       └── otu_table_L6.tsv
└── sample_list.txt
```
<br><br>
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
<br><br>
## Possible functions of the script
You can consult the function details of the pipeline invoking the main function name:

```bash
singularity exec K16S.sif kraken16S
```
![image](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/523df061-48a8-4a38-9907-ffdfa5b682d9)














