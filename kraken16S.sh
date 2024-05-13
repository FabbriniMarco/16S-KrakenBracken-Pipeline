#!/bin/bash

RED='\033[0;31m'
NC='\033[0m' # No Color
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
WHITE='\033[0;37m'
GREY='\033[1;30m'
LBL='\033[1;36m'
YELLOW='\033[0;33m'

# CONFIGURATION
DB_PATH="/mnt/databases/kraken_db"

check_command() {
    if ! command -v "$1" &> /dev/null; then
        return 1
    fi
    return 0
}

draw_progress_bar() {
    # Calculate the percentage completed
    percent=$((progress * 100 / tot))
    # Calculate the number of '#' characters to represent the progress
    num_chars=$((progress * 20 / tot))
    # Draw the progress bar
    printf "[%-20s] %d%%\r" "$(< /dev/zero tr '\0' '#' | head -c $num_chars)" "$percent"
}

# Check if kraken2 is available
if ! check_command "kraken2"; then
    echo "Error: kraken2 is not installed or cannot be accessed."
    exit 1
fi

# Check if bracken is available
if ! check_command "bracken"; then
    echo "Error: bracken is not installed or cannot be accessed."
    exit 1
fi

# Check if fastp is available
if ! check_command "fastp"; then
    echo "Error: fastp is not installed or cannot be accessed."
    exit 1
fi

# Function to log and display command output
log_and_display() {
    message="$1"
    # Log to the file without ANSI escape codes
    echo -e "$message" | sed -r 's/\x1B\[[0-9;]*[A-Za-z]//g' >> "$LOG_FILE"
    # Display to the screen with ANSI escape codes
    echo -e "$message"
}

# Get input data
if [ $# -eq 0 ]
then
    echo -e "\n${GREY}This script's path:${NC} $0

    ${BLUE}-i/--input${NC} : folder containing the .fastq files (paired- or/and single-end sequences) - ${RED}NO ABSOLUTE PATH${NC}
    ${BLUE}-p/--n_threads${NC} : number of processors used to run Kraken2 [${PURPLE}default: 10${NC}]
    ${BLUE}-P/--pigz_threads${NC} : number of processors used to run PIGZ for gzipping steps [${PURPLE}default: 10${NC}]
    ${BLUE}-D/--database${NC} : Database used: SILVA 138_1 SSURef NR99 [${PURPLE}default: silvaNR99${NC}]
    ${BLUE}-f/--filter_genome${NC} : filter the raw reads to remove host or contaminant genome using BMTagger with KneadData. Only for paired-end reads.
                Set to TRUE in order to use the database binded to /mnt/databases/bmtagger or add a local FULL PATH to a BMTagger indexed database [${PURPLE}default: FALSE${NC}]
    ${BLUE}-c/--check_fastp${NC} : [${GREEN}'TRUE'${NC} or ${GREEN}'FALSE'${NC}] specify whether to perform or not fastp error correction on RAW reads [${PURPLE}default: TRUE${NC}]
    ${BLUE}-r/--read_length_paired${NC} : [${GREEN}'250',  '500'${NC}] Set the read length number to select the proper version of Bracken database for the paired-end reads.
                [${PURPLE}default: 250${NC}]
    ${BLUE}-R/--read_length_single${NC} : [${GREEN}'250',  '500'${NC}]Set the read length number to select the proper version of Bracken database for the single-end reads.
                [${PURPLE}default: 500${NC}]
    ${BLUE}--no_diversity${NC} : If this flag is present, the script will only generate otu_tables and will not perform diversity analyses (only sample profiling)
    ${BLUE}--keeptemp${NC} : If this flag is present, the script will preserve intermediate files

    ${RED}
    # Single-end sequences MUST be in the ${samplename}_1.fastq format
    # There must always be at least a single paired-end sample to run the analysis
    ${NC}

    ${WHITE}Example usage:${NC}
    ${GREY} Full command: ${NC}
    ${LBL}kraken16S -i rawseqs -p 10 -P 10 -D silvaNR99 -f FALSE -c TRUE -r 250 -R 500${NC}
    ${GREY} Short command: ${NC}
    ${LBL}kraken16S -i rawseqs${NC}

    "
    exit 1
fi

# -----------------------
INPUT_FOLDER='';PROC='10';THREADS='10';DATABASE='silvaNR99';USEFASTP="TRUE";GENOMEFILTER="FALSE";RLPAIRED=250;RLSINGLE=500;CALCDIV=true;CLEAN=true;

OPTS=$(getopt -o i:p:P:D:f:c:r:R: -l input:,n_threads:,pigz_threads:,database:,filter_genome:,check_fastp:,read_length_paired:,read_length_single:,no_diversity,keeptemp, -- "$@")
eval set -- "$OPTS"

while [ $# -gt 0 ] ; do
    case "$1" in
    -i|--input) INPUT_FOLDER=$2; shift;;
    -p|--n_threads) PROC=$2; shift;;
    -P|--pigz_threads) THREADS=$2; shift;;
    -D|--database) DATABASE=$2; shift;;
    -f|--filter_genome) GENOMEFILTER=$2; shift;;
    -c|--check_fastp) USEFASTP=$2; shift;;
    -r|--read_length_paired) RLPAIRED=$2; shift;;
    -R|--read_length_single) RLSINGLE=$2; shift;;
    --no_diversity) CALCDIV=false;;
    --keeptemp) CLEAN=false ;; 
    (--) shift; break;;
    esac
    shift
done


# Check if input folder exists
if [ ! -d "$INPUT_FOLDER" ]; then
    echo -e "${RED}Error: Input folder does not exist.${NC}"
    exit 1
fi

start_time=`date +%s`
KVERSION=$(kraken2 --version | head -1)
BVERSION=$(bracken -v)

mkdir AnalysisKraken16_$(date +%F)_${DATABASE}
ls $INPUT_FOLDER | cut -d_ -f 1 | sed -n 'p;n' > AnalysisKraken16_$(date +%F)_${DATABASE}/sample_list.txt
SAMPLES_PROCESSED=$(wc -l AnalysisKraken16_$(date +%F)_${DATABASE}/sample_list.txt | cut -f1 -d " ");
cd AnalysisKraken16_$(date +%F)_${DATABASE}

# Specify the log file path
LOG_FILE="pipeline_log.txt"

if [ -n "$(LC_ALL=C type -t pigz)" ] && [ "$(LC_ALL=C type -t pigz)" = file ]; then 
    log_and_display "\n\t\t${GREY}pigz detected - proceeding with the script${NC}"; 
else 
    log_and_display -e "${RED}pigz is NOT installed${NC}.\nYou can download it from https://zlib.net/pigz/ - untar, enter folder and type 'make' " ; 
    exit 127 ; 
fi

if ls ../"$INPUT_FOLDER"/*.gz 1> /dev/null 2>&1; then
    log_and_display "${GREY}\t\t.gz extension detected${NC}"
    log_and_display "${GREY}\t\tDecompressing reads using pigz${NC}"
    pigz -p $THREADS -d ../$INPUT_FOLDER/*.gz
    log_and_display "${GREY}\t\tDecompression completed${NC}"
else if ls ../"$INPUT_FOLDER"/*.fastq 1> /dev/null 2>&1; then
    log_and_display "${GREY}\t\tAll files in the input folder are uncompressed .fastqs${GREEN}"
    log_and_display "${GREEN}\t\tProceeding${GREEN}"
fi ;
fi

# Diverging path for paired-end and single-end reads
# Read each sample from sample_list.txt
while IFS= read -r sample; do
    # Check if both _1.fastq and _2.fastq files exist for the sample
    if  { [ -e "../$INPUT_FOLDER/$sample"_1.fastq ] && [ -e "../$INPUT_FOLDER/$sample"_2.fastq ] ; } || { [ -e "../$INPUT_FOLDER/$sample"_*R1*.fastq ] && [ -e "../$INPUT_FOLDER/$sample"_*R2*.fastq ] ; } ; then
        #echo "Sample $sample has two files"
        echo ${sample} >> paired_sample_list.txt
    # Check if only one .fastq file exists for the sample
    elif [ -e "../$INPUT_FOLDER/$sample"_1.fastq ]; then
        #echo "Sample $sample has one file: ${sample}.fastq"
        echo ${sample} >> single_sample_list.txt
    else
        #echo "Sample $sample does not have the expected files";
        log_and_display "Sample $sample does not have the expected files."
        while true; do
            log_and_display "There might be an error in the sequence names or in the script. Please check!" ;
            read -p "Do you wish to proceed anyway? [Yy/Nn/QUIT]  " ANSWER
            sleep 1 ; echo ;
            case $ANSWER in
                [Yy]* ) log_and_display "\n\t\t${YELLOW}Proceeding with the script anyway${NC}"; sleep 2 ; break;;
                [Nn]* ) ;;
                [QUIT]* ) exit ;;
                * ) log_and_display "Please digit Y/y or N/n or digit QUIT to halt the script";;
            esac
        done
    fi
done < sample_list.txt

mkdir seqs

ls ../$INPUT_FOLDER/ | if grep -q  "_L001_R[0-2]_001.fastq"; then 
    for x in $(cat paired_sample_list.txt); do
    cp $(pwd)/../$INPUT_FOLDER/${x}_*R1*.fastq seqs/${x}_1.fastq ;
    cp $(pwd)/../$INPUT_FOLDER/${x}_*R2*.fastq seqs/${x}_2.fastq ;
    done;
else 
    for x in $(cat paired_sample_list.txt); do
    cp $(pwd)/../$INPUT_FOLDER/${x}_*1.fastq seqs/${x}_1.fastq ;
    cp $(pwd)/../$INPUT_FOLDER/${x}_*2.fastq seqs/${x}_2.fastq ;
    done  ;
fi ;

if  [ -e "single_sample_list.txt" ]; then
    for x in $(cat single_sample_list.txt); do
    cp $(pwd)/../$INPUT_FOLDER/${x}_1.fastq seqs/${x}_1.fastq ;
    done
fi

for x in $(cat sample_list.txt); do
wc -l seqs/${x}_1.fastq >> reads_counts.txt ;
done

#### FASTP error correction
if [ "$USEFASTP" = "TRUE" ]; then
    log_and_display "\n\t\t${GREY}Performing read correction with fastp${NC}";
    mkdir seqs_corrected
    samplen=1
    progress=0
    tot=$(wc -l < paired_sample_list.txt)
    for x in $(cat paired_sample_list.txt); do
        fastp -w ${PROC} -i seqs/${x}_1.fastq -I seqs/${x}_2.fastq -o seqs_corrected/${x}_1.fastq -O seqs_corrected/${x}_2.fastq 2> /dev/null ;
        ((samplen++)) ;
        ((progress++));
        draw_progress_bar ;
    done
    if  [ -e "single_sample_list.txt" ]; then
        samplen=1
        progress=0
        tot=$(wc -l < single_sample_list.txt)
        for x in $(cat single_sample_list.txt); do
        fastp -w ${PROC} -i seqs/${x}_1.fastq -o seqs_corrected/${x}_1.fastq 2> /dev/null ;
        ((samplen++)) ;
        ((progress++));
        draw_progress_bar ;
        done
    fi
    mv seqs seqs_uncorrected
    mv seqs_corrected seqs
    rm -fr fastp.json fastp.html
    log_and_display "\n\t\t${GREEN}Read correction with fastp completed${NC}";
fi


### Genome filtering if requested
if [ "$GENOMEFILTER" = "FALSE" ]; then
    log_and_display "${GREY}\n\t\tNo genome filtering requested${NC}";
else if [ "$GENOMEFILTER" = "TRUE" ]; then
    log_and_display "\n\t\t${GREY}Performing Genome filtering with binded database${NC}";
    mkdir seqs_filtered
    GENOMEPATH="/mnt/databases/bmtagger/"
    # Filter out human reads using the binded database
    for x in $(cat paired_sample_list.txt); do
    ((i=i%${PROC})); ((i++==0)) && wait
        kneaddata --input seqs/${x}_1.fastq --input seqs/${x}_2.fastq \
        --run-bmtagger \
        --reference-db ${GENOMEPATH} \
        --output seqs_filtered/ \
        --threads 1 --processes 1 \
        --max-memory 1g \
        --bypass-trim \
        --remove-intermediate-output &
    done
    wait
    # Fix headers
    for x in $(cat paired_sample_list.txt); do 
        sed -i "s/\\#0\\/1//g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ; 
        sed -i "s/:N/ 1:N/g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ;
        sed -i "s/\\#0\\/2//g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
        sed -i "s/:N/ 2:N/g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
        mv seqs_filtered/${x}_1_kneaddata_paired_1.fastq seqs_filtered/${x}_1.fastq ;
        mv seqs_filtered/${x}_1_kneaddata_paired_2.fastq seqs_filtered/${x}_2.fastq ;
    done
    for x in $(cat paired_sample_list.txt); do
        wc -l seqs_filtered/${x}_1.fastq >> reads_counts_FILTERED.txt ;
    done
    if  [ -e "single_sample_list.txt" ]; then
        for x in $(cat single_sample_list.txt); do
        ((i=i%${PROC})); ((i++==0)) && wait
            kneaddata --input seqs/${x}_1.fastq \
            --run-bmtagger \
            --reference-db ${GENOMEPATH} \
            --output seqs_filtered/ \
            --threads 1 --processes 1 \
            --max-memory 1g \
            --bypass-trim \
            --remove-intermediate-output &
        done
        wait
        # Fix headers
        for x in $(cat single_sample_list.txt); do 
            sed -i "s/\\#0\\/1//g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ; 
            sed -i "s/:N/ 1:N/g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ;
            sed -i "s/\\#0\\/2//g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
            sed -i "s/:N/ 2:N/g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
            mv seqs_filtered/${x}_1_kneaddata_paired_1.fastq seqs_filtered/${x}_1.fastq ;
            mv seqs_filtered/${x}_1_kneaddata_paired_2.fastq seqs_filtered/${x}_2.fastq ;
        done
        for x in $(cat single_sample_list.txt); do
            wc -l seqs_filtered/${x}_1.fastq >> reads_counts_FILTERED.txt ;
        done
    fi
    mv seqs seqs_unfiltered
    mv seqs_filtered seqs
    log_and_display "\n\t\t${GREEN}Genome filtering completed with genome: ${GENOMEPATH}${NC}";
else if echo $GENOMEFILTER | grep -q "/" ; then
    mkdir seqs_filtered
    #Check if directory exists
    if [ -d "$GENOMEFILTER" ]; then
        log_and_display "\n\t\t${GREEN}Custom directory provided found in: ${GENOMEFILTER}${NC}"
        # Filter out human reads using the provided custom database
        for x in $(cat paired_sample_list.txt); do
        ((i=i%${PROC})); ((i++==0)) && wait
            kneaddata --input seqs/${x}_1.fastq --input seqs/${x}_2.fastq \
            --run-bmtagger \
            --reference-db ${GENOMEFILTER} \
            --output seqs_filtered/ \
            --threads 1 --processes 1 \
            --max-memory 1g \
            --bypass-trim \
            --remove-intermediate-output &
        done
        wait
        # Fix headers
        for x in $(cat paired_sample_list.txt); do 
            sed -i "s/\\#0\\/1//g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ; 
            sed -i "s/:N/ 1:N/g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ;
            sed -i "s/\\#0\\/2//g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
            sed -i "s/:N/ 2:N/g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
            mv seqs_filtered/${x}_1_kneaddata_paired_1.fastq seqs_filtered/${x}_1.fastq ;
            mv seqs_filtered/${x}_1_kneaddata_paired_2.fastq seqs_filtered/${x}_2.fastq ;
        done
        for x in $(cat paired_sample_list.txt); do
            wc -l seqs_filtered/${x}_1.fastq >> reads_counts_FILTERED.txt ;
        done
        if  [ -e "single_sample_list.txt" ]; then
            for x in $(cat single_sample_list.txt); do
            ((i=i%${PROC})); ((i++==0)) && wait
                kneaddata --input seqs/${x}_1.fastq \
                --run-bmtagger \
                --reference-db ${GENOMEPATH} \
                --output seqs_filtered/ \
                --threads 1 --processes 1 \
                --max-memory 1g \
                --bypass-trim \
                --remove-intermediate-output &
            done
            wait
            # Fix headers
            for x in $(cat single_sample_list.txt); do 
                sed -i "s/\\#0\\/1//g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ; 
                sed -i "s/:N/ 1:N/g" seqs_filtered/${x}_1_kneaddata_paired_1.fastq ;
                sed -i "s/\\#0\\/2//g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
                sed -i "s/:N/ 2:N/g" seqs_filtered/${x}_1_kneaddata_paired_2.fastq ;
                mv seqs_filtered/${x}_1_kneaddata_paired_1.fastq seqs_filtered/${x}_1.fastq ;
                mv seqs_filtered/${x}_1_kneaddata_paired_2.fastq seqs_filtered/${x}_2.fastq ;
            done
            for x in $(cat single_sample_list.txt); do
                wc -l seqs_filtered/${x}_1.fastq >> reads_counts_FILTERED.txt ;
            done
        fi
        mv seqs seqs_unfiltered
        mv seqs_filtered seqs
        log_and_display "\n\t\t${GREEN}Genome filtering completed with genome: ${GENOMEFILTER}${NC}";
    fi ;
else log_and_display "${RED}Wrong parameter specified for '-f/--filter_genome'; remember to use ${BLUE}TRUE${RED} or ${BLUE}FALSE${RED} or a ${BLUE}FULL PATH${RED} to an indexed BMTagger database${NC}"; exit 2; 
fi;
fi;
fi

### Kraken2 mapping
mkdir krakenized
log_and_display "\n\t\t${GREY}Performing Kraken2 read mapping against ${DATABASE}${NC}";

if  [ -e "single_sample_list.txt" ]; then
    samplen=1
    progress=0
    tot=$(wc -l < single_sample_list.txt)
    for i in $(cat single_sample_list.txt); do 
        kraken2 --db ${DB_PATH}/${DATABASE} \
        --threads ${PROC} \
        --report krakenized/${i}.kreport seqs/${i}_1.fastq > krakenized/${i}.kraken 2> /dev/null ;
        ((samplen++)) ;
        ((progress++));
        draw_progress_bar ;
    done
fi

samplen=1
progress=0
tot=$(wc -l < paired_sample_list.txt)
for i in $(cat paired_sample_list.txt); do 
    kraken2 --db ${DB_PATH}/${DATABASE} \
    --threads ${PROC} \
    --report krakenized/${i}.kreport \
    --paired seqs/${i}_1.fastq seqs/${i}_2.fastq > krakenized/${i}.kraken 2> /dev/null ;
    ((samplen++)) ;
    ((progress++));
    draw_progress_bar ;
done

log_and_display "\n\t\t${GREEN}Kraken2 mapping completed${NC}";


## Bracken2 mapping
mkdir brackenized

log_and_display "\n\t\t${GREY}Running Bracken with ${RLPAIRED}bp and ${RLSINGLE}bp parameters${NC}";

samplen=1
progress=0
tot=$(wc -l < paired_sample_list.txt)
for i in $(cat paired_sample_list.txt); do 
    bracken -d ${DB_PATH}/${DATABASE} -i krakenized/${i}.kreport -o brackenized/${i}_phyla.bracken -r ${RLPAIRED} -l P 1> /dev/null 2> /dev/null ;
    bracken -d ${DB_PATH}/${DATABASE} -i krakenized/${i}.kreport -o brackenized/${i}_families.bracken -r ${RLPAIRED} -l F 1> /dev/null 2> /dev/null ;
    bracken -d ${DB_PATH}/${DATABASE} -i krakenized/${i}.kreport -o brackenized/${i}_genera.bracken -r ${RLPAIRED} -l G 1> /dev/null 2> /dev/null ;
    ((samplen++)) ;
    ((progress++));
    draw_progress_bar ;
done

if  [ -e "single_sample_list.txt" ]; then
    samplen=1
    progress=0
    tot=$(wc -l < single_sample_list.txt)
    for i in $(cat single_sample_list.txt); do 
        bracken -d ${DB_PATH}/${DATABASE} -i krakenized/${i}.kreport -o brackenized/${i}_phyla.bracken -r ${RLSINGLE} -l P 1> /dev/null ;
        bracken -d ${DB_PATH}/${DATABASE} -i krakenized/${i}.kreport -o brackenized/${i}_families.bracken -r ${RLSINGLE} -l F 1> /dev/null ;
        bracken -d ${DB_PATH}/${DATABASE} -i krakenized/${i}.kreport -o brackenized/${i}_genera.bracken -r ${RLSINGLE} -l G 1> /dev/null ;
        ((samplen++)) ;
        ((progress++));
        draw_progress_bar ;
    done
fi 

log_and_display "\n\t\t${GREEN}Bracken completed. Genearting otu tables${NC}";

# Gather the reports in an otu table
samplen=1
progress=0
tot=$(wc -l < sample_list.txt)
for i in $(cat sample_list.txt); do
    awk -v OFS=$'\t' -v sample=${i} '{print sample, $0}' <( sed '1d' brackenized/${i}_phyla.bracken) >> otu_table_final.tsv
    awk -v OFS=$'\t' -v sample=${i} '{print sample, $0}' <( sed '1d' brackenized/${i}_families.bracken) >> otu_table_final.tsv
    awk -v OFS=$'\t' -v sample=${i} '{print sample, $0}' <( sed '1d' brackenized/${i}_genera.bracken) >> otu_table_final.tsv
    ((samplen++)) ;
    ((progress++))
    draw_progress_bar
done && sed -i '1 i\ID\ttaxa\ttaxid\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads' otu_table_final.tsv

log_and_display "\n\t\t${GREY}Formatting otu tables${NC}"
CPU_process_otutab_relabb.py
CPU_process_otutab.py
log_and_display "\n\t\t${GREEN}Otu tables generated${NC}"

mkdir RESULTS
mv otu_table_final.tsv RESULTS/.
mv final_otu_table_matrix_relabb.csv RESULTS/.
mv final_otu_table_matrix.csv RESULTS/.
mv taxa_summary RESULTS/.

if $CALCDIV; then
    log_and_display "\n\t\t${GREY}Computing alpha and beta diversity metrics${NC}"
    cd RESULTS
    compute_diversity_kraken.r 
    cd ..
    log_and_display "\n\t\t${GREEN}Diversity analyses completed${NC}"
fi

if $CLEAN; then
    rm -fr seqs*
    rm -fr brackenized
    rm -fr krakenized
fi

log_and_display "\n\t\t${GREY}Compressing reads...${NC}"
pigz -p 9 ../${INPUT_FOLDER}/*.fastq;

if $CLEAN; then
    find . -type f -name '*.fastq' -exec pigz -p 9 {} +
fi

log_and_display "\n\t\t${GREEN}Compression complete.${NC}"

end_time=`date +%s`
elapsed=$(( (end_time-start_time)/60 ))

##-----------

log_and_display "${GREEN}Done${NC}\n${YELLOW}---> Check your data${NC}"

log_and_display " \nParameter recap:
    ~ samples processed: ${SAMPLES_PROCESSED}
    ~ contaminant genome filtered: ${GENOMEFILTER}
    ~ database: ${DATABASE}
    ~ Kraken2 threads: ${PROC}
    ~ pigz threads: ${THREADS}
    ~ fastp read correction : ${USEFASTP}
    ~ Read length for Bracken (paired/single): ${RLPAIRED}/${RLSINGLE} bp
    ~ Kraken2 version: ${KVERSION}
    ~ Bracken version: ${BVERSION}
    
    Time required: ${elapsed} minutes" > RESULTS/logs.log 

log_and_display " \nParameter recap:
    ~ samples processed: ${SAMPLES_PROCESSED}
    ~ contaminant genome filtered: ${GENOMEFILTER}
    ~ database: ${DATABASE}
    ~ Kraken2 threads: ${PROC}
    ~ pigz threads: ${THREADS}
    ~ fastp read correction : ${USEFASTP}
    ~ Read length for Bracken (paired/single): ${RLPAIRED}/${RLSINGLE} bp
    ~ Kraken2 version: ${KVERSION}
    ~ Bracken version: ${BVERSION}

    Time required: ${elapsed} minutes"
