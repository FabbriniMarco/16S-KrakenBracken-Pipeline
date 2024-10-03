#!/usr/bin/env nextflow

// Default parameters
// params.output_folder = "./AnalysisKraken_${currentDate}"          // Output folder
// params.fastq_folder = "seqs"                                      // Folder with FASTQ files
// params.genome_path = "/mnt/luks/databases/hg38_bowtie2"           // Path to host filter database
// params.kraken_db_path = "/mnt/databases/kraken_db/silvaNR99"      // Path to kraken2 silvaNR99 database
// params.readlength_paired = 250                                    // Read length for bracken database for paired-end reads
// params.readlength_single = 500                                    // Read length for bracken database for single-end reads
// params.use_fastp = true                                           // Use fastp to clean and clip reads
// params.use_kneaddata = true                                       // Perform host filtering using KneadData implementing bowtie2
// params.help = false                                               // Print the help

// Display help message if params.help is true
if (params.help) {
    log.info "\u001B[1;34mUsage:\u001B[0m nextflow run kraken16S.nf [options]\n" +
    """
    \u001B[1;32mOptions:\u001B[0m
      \u001B[1;33m--fastq_folder\u001B[0m        Path to the folder containing input FASTQ files [\u001B[1;35mdefault: './seqs'\u001B[0m]
      \u001B[1;33m--output_folder\u001B[0m       Path to the folder for storing results [\u001B[1;35mdefault: './AnalysisKraken_\${currentDate}'\u001B[0m]
      \u001B[1;33m--genome_path\u001B[0m         Path to the genome reference database for KneadData [\u001B[1;35mdefault: 'none'\u001B[0m, MUST be provided if use_kneaddata is set to true]
      \u001B[1;33m--kraken_db_path\u001B[0m      Path to the Kraken2 database [\u001B[1;35mdefault: '/mnt/databases/kraken_db/silvaNR99'\u001B[0m inside the singularity container]
      \u001B[1;33m--readlength_paired\u001B[0m   Read length for paired-end reads (250/500) [\u001B[1;35mdefault: 250\u001B[0m]
      \u001B[1;33m--readlength_single\u001B[0m   Read length for single-end reads (250/500) [\u001B[1;35mdefault: 500\u001B[0m]
      \u001B[1;33m--use_fastp\u001B[0m           Enable or disable fastp (true/false) [\u001B[1;35mdefault: true\u001B[0m]
      \u001B[1;33m--use_kneaddata\u001B[0m       Enable or disable KneadData (true/false) [\u001B[1;35mdefault: true\u001B[0m]
      \u001B[1;33m--help\u001B[0m                Show this help message

    \u001B[1;32mExample fast run (use max processor available and skip host genome filtering):\u001B[0m
      \u001B[1;36mnextflow run main.nf --fastq_folder my_raw_seqs --use_kneaddata false\u001B[0m

    \u001B[1;32mExample with fine host genome filtering and parallelization control:\u001B[0m
      \u001B[1;36mnextflow run main.nf -process.cpus 40 -process.memory '128 GB' -process.maxForks 30 --fastq_folder raw_seqs --genome_path HOST_BOWTIE2_DB_PATH\u001B[0m

    \u001B[1;32mTo clean working file once the analysis is completed:\u001B[0m
      \u001B[1;36mnextflow clean -f -q\u001B[0m

    \u001B[1;32mIf you don't need them you can also delete the folders:\u001B[0m
      \u001B[1;30m'bracken', 'fastp', 'kraken2' and 'host_filtering' inside the ./AnalysisKraken_\${currentDate} folder\u001B[0m
    """
    exit 0
}


// Channel for paired-end files (compressed and uncompressed)
Channel
    .fromFilePairs("${params.fastq_folder}/*_{1,2}.fastq*")
    .set { paired_end_files }

// Channel for single-end files (compressed and uncompressed)
Channel
    .fromPath("${params.fastq_folder}/*.fastq*")
    .filter { !it.getName().matches(/.*_[12]\.fastq.gz$/) }
    .map { file -> tuple(file.getSimpleName().replaceAll(/\.fastq(\.gz)?$/, ''), file) }
    .set { single_end_files }

// Merge paired and single-end files
all_fastq_files = paired_end_files.mix(single_end_files)

// Filter compressed (.fastq.gz) and uncompressed (.fastq) files
compressed_fastq_files = all_fastq_files.filter { it -> it[1] instanceof List ? it[1].every { f -> f.name.endsWith('.fastq.gz') } : it[1].name.endsWith('.fastq.gz') }
uncompressed_fastq_files = all_fastq_files.filter { it -> it[1] instanceof List ? it[1].every { f -> f.name.endsWith('.fastq') } : it[1].name.endsWith('.fastq') }


// Decompress only compressed (.fastq.gz) files
process compress_reads {
    cleanup 'onSuccess'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}*.fastq.gz")

    stageInMode 'copy'

    script:
    """
    pigz -p 1 ${reads}
    """
}

// Define the process for running fastp
process fastp {
    //cleanup 'onSuccess'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}/*.fastq.gz")

    publishDir "${params.output_folder}/fastp", mode: 'copy'
    
    script:
    """
    mkdir -p ${id}
    # Detect if it's paired-end or single-end based on the number of elements in reads
    if [ -f "${reads[1]}" ]; then
        # Paired-end data
        fastp -i ${reads[0]} -I ${reads[1]} -o ${id}/${id}_1.fastq.gz -O ${id}/${id}_2.fastq.gz
    else
        # Single-end data
        fastp -i ${reads[0]} -o ${id}/${id}.fastq.gz
    fi
    """
}


// Define the process for running kneaddata for host filtering
process kneaddata {
    //cleanup 'onSuccess'

    publishDir "${params.output_folder}/host_filtering", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}/*.fastq.gz")

    //stageInMode 'copy'

    script:
    """
    mkdir -p ${id}
    # Detect if it's paired-end or single-end based on the number of elements in reads
    if [ -f "${reads[1]}" ]; then
        # Paired-end data
        kneaddata --input ${reads[0]} --input ${reads[1]} \\
        --reference-db ${params.genome_path} \\
        --output . \\
        --threads 1 --processes 1 \\
        --max-memory 3g \\
        --bypass-trim \\
        --remove-intermediate-output

        # Adjust headers and rename files
        rm -fr reads *_hg38_* *unmatched*
        mv ${id}_*_paired_1.fastq ${id}/${id}_1.fastq
        mv ${id}_*_paired_2.fastq ${id}/${id}_2.fastq
        sed -i "s/\\#0\\/1//g" ${id}/${id}_1.fastq
        sed -i "s/:N/ 1:N/g" ${id}/${id}_1.fastq
        sed -i "s/\\#0\\/2//g" ${id}/${id}_2.fastq
        sed -i "s/:N/ 2:N/g" ${id}/${id}_2.fastq
        pigz -p 3 ${id}/*.fastq
    else
        # Single-end data
        kneaddata --input ${reads[0]} \\
        --reference-db ${params.genome_path} \\
        --output . \\
        --threads 1 --processes 1 \\
        --max-memory 3g \\
        --bypass-trim \\
        --remove-intermediate-output

        # Adjust headers
        rm -fr reads *_hg38_*
        mv ${id}_kneaddata.fastq ${id}/${id}.fastq
        sed -i "s/\\#0\\/1//g" ${id}/${id}.fastq
        sed -i "s/:N/ 1:N/g" ${id}/${id}.fastq
        pigz -p 3 ${id}/${id}.fastq
    fi
    """
}


// Define the process for running kraken2 mapping
process kraken2 {
    cleanup 'onSuccess'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}/${id}_*.kreport")

    publishDir "${params.output_folder}/kraken2", mode: 'copy'

    //stageInMode 'copy'
    
    script:
    """
    mkdir -p ${id}
    # Detect if it's paired-end or single-end based on the number of elements in reads
    if [ -f "${reads[1]}" ]; then
        # Paired-end data
        kraken2 --db ${params.kraken_db_path} --gzip-compressed --threads 3 --report ${id}/${id}_paired-end.kreport --paired ${reads[0]} ${reads[1]} 1> /dev/null 2> /dev/null
    else
        # Single-end data
        kraken2 --db ${params.kraken_db_path} --gzip-compressed --threads 3 --report ${id}/${id}_single-end.kreport ${reads} 1> /dev/null 2> /dev/null
    fi
    """
}


process bracken {
    cleanup 'onSuccess'

    input:
    tuple val(id), path(kreportFile)

    output:
    tuple val(id), path("${id}/${id}_*.bracken")

    publishDir "${params.output_folder}/bracken", mode: 'copy'

    stageInMode 'copy'

    script:
    """
    mkdir -p ${id}
    # Detect if it's paired-end or single-end based on the number of elements in reads
    if [[ "${kreportFile}" == *"_paired-end.kreport" ]]; then
        # Paired-end data
        bracken -d ${params.kraken_db_path} -i ${id}_paired-end.kreport -o ${id}/${id}_phyla.bracken -r ${params.readlength_paired} -l P
        bracken -d ${params.kraken_db_path} -i ${id}_paired-end.kreport -o ${id}/${id}_families.bracken -r ${params.readlength_paired} -l F
        bracken -d ${params.kraken_db_path} -i ${id}_paired-end.kreport -o ${id}/${id}_genera.bracken -r ${params.readlength_paired} -l G
    else
        # Single-end data
        bracken -d ${params.kraken_db_path} -i ${id}_single-end.kreport -o ${id}/${id}_phyla.bracken -r ${params.readlength_single} -l P
        bracken -d ${params.kraken_db_path} -i ${id}_single-end.kreport -o ${id}/${id}_families.bracken -r ${params.readlength_single} -l F
        bracken -d ${params.kraken_db_path} -i ${id}_single-end.kreport -o ${id}/${id}_genera.bracken -r ${params.readlength_single} -l G
    fi
    """
}

process generate_otu_tables {
    input:
    tuple val(id), path(kreportFiles)

    output:
    path "${id}_otu_table.tsv"

    script:
    firstFile = kreportFiles[0]
    secondFile = kreportFiles[1]
    thirdFile = kreportFiles[2]

    """
    # Append families, genera, and phyla to the ${id}_otu_table.tsv
    awk -v OFS='\t' -v sample=${id} '{print sample, \$0}' <(sed '1d' ${thirdFile}) >> ${id}_otu_table.tsv
    awk -v OFS='\t' -v sample=${id} '{print sample, \$0}' <(sed '1d' ${firstFile}) >> ${id}_otu_table.tsv
    awk -v OFS='\t' -v sample=${id} '{print sample, \$0}' <(sed '1d' ${secondFile}) >> ${id}_otu_table.tsv
    """
}


process combine_otu_tables {
    input:
    path otu_tables

    output:
    path "otu_table_final.tsv"

    publishDir "${params.output_folder}/RESULTS", mode: 'copy'

    script:
    """
    cat ${otu_tables} > otu_table_final.tsv
    sed -i "1 i\\ID\ttaxa\ttaxid\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads" otu_table_final.tsv
    """
}


process generate_matrix_otutables {
    input:
    path combined_otu_table

    output:
    tuple path("*.csv"), path("taxa_summary/")

    publishDir "${params.output_folder}/RESULTS", mode: 'copy'

    script:
    """
    CPU_process_otutab.py
    CPU_process_otutab_relabb.py
    """
}


process calculate_diversity {
    input:
    tuple path(final_otu_table), path(final_otu_relabb)

    output:
    tuple path("beta_diversity/"), path("alpha.xlsx")

    publishDir "${params.output_folder}/RESULTS", mode: 'copy'

    script:
    """
    compute_diversity_kraken.r
    """
}


// Define the workflow
workflow {
    // Compress eventually decompressed .fastq files
    compressed_files = compress_reads(uncompressed_fastq_files)

    // Merge the freshly compressed files with the already compressed ones
    all_compressed_files = compressed_files.mix(compressed_fastq_files)

    // Process reads with fastp if requested
    if (params.use_fastp) {
        fastp_out = fastp(all_compressed_files)
    } else {
        fastp_out = all_compressed_files
    }

    // Filter host genome with KneadData if requested
    if (params.use_kneaddata) {
        kneaddata_out = kneaddata(fastp_out)
    } else {
        kneaddata_out = fastp_out
    }

    // Align reads on Silva SSURef 138.1 NR99 with Kraken2
    kraken2_out = kraken2(kneaddata_out)

    // Collect compositional details with Bracken
    bracken_out = bracken(kraken2_out)

    // Generate a tidy otu table
    otutab = generate_otu_tables(bracken_out)
    otutab_list = otutab.collect()

    combined_otu_table = combine_otu_tables(otutab_list)
    matrix_otutables = generate_matrix_otutables(combined_otu_table)
    diversity = calculate_diversity(matrix_otutables)
}
// 03.10.2024 - Fabbrini Marco - fabbrinimarcoo.mf@gmail.com