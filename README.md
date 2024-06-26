# 16S-KrakenBracken-Pipeline

This repository hosts a containerized pipeline tailored for 16S rRNA amplicon sequencing analyses, leveraging [Kraken2](https://github.com/DerrickWood/kraken2)/[Bracken](https://github.com/jenniferlu717/Bracken) to generate abundance tables (OTU tables) and compute alpha and beta diversity results.

The pipeline is designed to process FASTQ 16S amplicon files, preferably in paired-end format. While it can handle mixed paired-end and single-end data, it does not support solely single-end data. Upon input, sequences undergo quality filtering using [fastp](https://github.com/OpenGene/fastp), ensuring high-quality data for downstream analysis. Additionally, the pipeline offers the capability to filter sequences for a reference genome using [KneadData](https://github.com/biobakery/kneaddata), providing flexibility for various research needs.

Processed sequences are then mapped onto the [SILVA 138_1 SSURef database](https://www.arb-silva.de/no_cache/download/archive/release_138_1/) using Kraken2, with taxonomic assignments obtained through Bracken. The resulting taxonomic outputs are further processed to generate counts and relative percentages for OTU tables. Additionally, the pipeline calculates alpha diversity values using several metrics and produces beta diversity distance matrices, enriching the depth of analysis.
A schematic overview of the pipeline is shown in the chart below, created in Lucidchart ([Lucid.co](https://lucid.co/)):

![Screenshot 2024-05-15 120420](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/e91ae0d9-3113-4dd4-bf58-86926a7d9076)


## Table of contents
1. [Dependencies](#dependencies)
2. [Installation](#installation)
	- [Pre-built container (recommended)](#pre-built-container-recommended)
	- [Build your own container](#build-your-own-container)
3. [Usage instructions](#usage-instruction)
	- [Running inside the container](#running-inside-the-container)
	- [Pass a script to the container](#pass-a-script-to-the-container)
4. [Pipeline outputs](#pipeline-outputs)
5. [Generating a BMTagger-compatible reference genome](#generating-a-bmtagger-compatible-reference-genome)
6. [Parameters of the script](#parameters-of-the-script)
7. [Citation](#citation)

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
singularity pull --arch amd64 library://fabbrinimarco/16s-krakenbracken-pipeline/k16s:sha256.7dc4939abc3f620d47e06edeba9adea39d8975f2504ee88b74642a209910d917
# Check the container
singularity run k16s_v1.sif kraken16S 
```

<br>

### Build your own container
Clone this repository or download the build_container.def file and use it to build the Singularity container:
```bash
singularity build k16s_v1.sif build_container.def
```
Note that during the Kraken2 Silva database building process, 24 threads will be utilized.
Due to the unconsistent access to NCBI FTP services, the container build might exit while downloading NCBI BLAST+.
If you encounter this issue and still want to build your own container you can build it without BLAST+ as sandbox using the kraken16S_container_noBLAST+.def from this repository, install BLAST+ manually and then recreate the SIF image:

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

## Usage instruction
### Running inside the container
Executing commands within the container is straightforward. Assuming you have a folder named 'rawseqs' containing your raw 16S FASTQ sequences (in .fastq or .fastq.gz format) in the current path, run:

``` bash
singularity exec --bind $(pwd):/workdir k16s_v1.sif bash
cd /workdir
kraken16S -i rawseqs/
# Once the script is completed remember to close the container shell
exit
```

You can also perform reference genome filtering, such as with the human genome, by binding a local BMTagger-indexed genome to the container folder '/mnt/databases/bmtagger', as follows:

``` bash
# In this example I have a BMTagger-indexed database in the local path '/mnt/luks/databases/hg38'
singularity exec --bind /mnt/luks/databases/hg38:/mnt/databases/bmtagger --bind $(pwd):/workdir k16s_v1.sif bash
cd /workdir
kraken16S -i rawseqs/ -f TRUE
```
<br>

### Pass a script to the container
You can also pass the container a shell script to be executed to avoid opening a shell inside the container:

```bash
echo -e "#!/bin/bash
cd /workdir
kraken16S -i rawseqs/" >> kraken16Scall.sh

chmod +x kraken16Scall.sh

singularity exec --bind $(pwd):/workdir k16s_v1.sif kraken16Scall.sh
```

Or with genome filtering, as the example above:

```bash
echo -e "#!/bin/bash
cd /workdir
kraken16S -i rawseqs/ -f TRUE" >> kraken16Scall.sh

chmod +x kraken16Scall.sh

singularity exec --bind /mnt/luks/databases/hg38:/mnt/databases/bmtagger --bind $(pwd):/workdir k16s_v1.sif kraken16Scall.sh
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
singularity exec k16s_v1.sif kraken16S
```


![Screenshot 2024-05-13 175759](https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline/assets/83694005/478a84ea-c17f-4b21-a6b4-cab3a0d8f2dc)


For example, you can control the number of parallel threads used for the read mapping and read compress/decompress operations. You can decide wether or not to perform fastp read quality check and filtering, as well as eventual host genome filtering. If you deal with sequences with different length rather than 250bp for paired-end data and 500bp for single-end data you can also edit the read length Bracken settings choosing either 250 or 500bp options.
You can also decide to perform only taxonomic assignment and skip the alpha dn beta diversity analyses, or keep the temp files (i.e., the kraken2 and bracken outputs, as well as filtered reads).

<br>

## Citation

You can cite the use of this pipeline as such:

```diff
@Manual{,
  title = {16S-KrakenBracken-Pipeline: Containerized pipeline for 16S rRNA amplicon sequencing analyses using Kraken2/Bracken},
  author = {Marco Fabbrini},
  year = {2024},
  note = {Singularity image v1},
  url = {https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline.git},
}
```
Example: 
> Fabbrini, M. (2024) '16S-KrakenBracken-Pipeline - Containerized pipeline for 16S rRNA amplicon sequencing analyses using Kraken2/Bracken'. Available on Github: https://github.com/FabbriniMarco/16S-KrakenBracken-Pipeline














