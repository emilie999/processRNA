#!/bin/bash


# help message
display_help()
{
    echo -e "`cat ./docs/help.txt`"
}


# detailed documentation
display_detailed_doc()
{
    display_help
    echo -e "`cat ./docs/doc.txt`"
}


# initialize default variables
INPUT_DIR=0
OUTPUT_DIR=0
MERGING=0
GENOME=0
VE=$HOME


# get options
while getopts "hdi:o:e:mpg:c:" option; do
    case $option in
        # display help message
        h)  display_help
            exit 0;;
        # display detailed documentation
        d)  display_detailed_doc
            exit 0;;
        # input directory
        i)  INPUT_DIR=$OPTARG
            if [[ -d $INPUT_DIR ]]; then
                echo "Reading input from $INPUT_DIR"
            else
                echo "ERROR: please enter a valid directory for option -i"
                exit 2
            fi;;
        # output directory
        o)  OUTPUT_DIR=$OPTARG
            if [[ -d $OUTPUT_DIR ]]; then
                echo "Writing output to $OUTPUT_DIR"
            else
                echo "ERROR: please enter a valid directory for option -o"
                exit 2
            fi;;
        # virtual environment
        e)  VE=$OPTARG
            if [[ -d $VE ]]; then
                echo "Using the venv $VE"
            else
                echo "ERROR: please enter a valid directory for option -e"
                exit 2
            fi;;
        # merging option (lanes)
        m)  MERGING=1;;
        # paired end option (reads)
        p)  PAIRED=1;;
        # reference genome
        g)  GENOME=$OPTARG
            if (($GENOME == "mm10")); then
                GENDIR=/mnt/data/genomes/mm10_star_new/
            elif (($GENOME == "hg38")); then
                GENDIR=/mnt/data/genomes/hg38_star_2.7.9a/
            else
                echo "ERROR: please enter a valid reference genome for option -g"
                exit 2
            fi;;
        # comparisons
        g)  COMPARISONS=$OPTARG;; # INCOMPLETE !!!!!!!!!!
        # exit on invalid options
        \?) echo "ERROR: invalid option usage"
            display_help
            exit 1;;
    esac
done


# check that the required arguments are given
if ! [[ -d $INPUT_DIR ]]; then
    echo "ERROR: missing input directory"
    display_help
    exit 1
fi
if ! [[ -f ${INPUT_DIR}/metadata.txt ]]; then
    echo "ERROR: missing metadata file in input directory"
    display_help
    exit 1
fi
if ! [[ -d $OUTPUT_DIR ]]; then
    echo "ERROR: missing output directory"
    display_help
    exit 1
fi
if (($GENOME == 0)); then
    echo "ERROR: missing genome"
    display_help
    exit 1
fi


$METADATA=${INPUT_DIR}/metadata.txt


# make subdirectories for outputs
qc=${OUTPUT_DIR}/QC
mapped=${OUTPUT_DIR}/mapped
tagir=${OUTPUT_DIR}/tagdir
counts=${OUTPUT_DIR}/counts
results=${OUTPUT_DIR}/results
mkdir $qc $mapped $tagdir $counts $results


# generate FASTQC results
fastqc -t 6 -q *.fastq.gz -o $qc


# make a list of unique sample names
samples=($(tail -n +2 "$METADATA" | awk -F '\t' '{print $1}'))
filenames=($(tail -n +2 "$METADATA" | awk -F '\t' '{print $2}'))
num=${#samples[@]}


# map to ref genome
for f in ${INPUT_DIR}/*.fastq.gz; do
    name=$(echo ${f##*/})
    subname=$(echo $name | cut -f1 -d.)
    STAR --genomeDir $GENDIR --runThreadN 16 --readFilesType Fastx --readFilesIn $f --readFilesCommand gunzip -c --outFileNamePrefix ${mapped}/${subname} --outSAMtype BAM SortedByCoordinate
done


# merging...
if (($MERGING==1)); then
    merged=${OUTPUT_DIR}/merged
    mkdir $merged
    for sample in "${samples[@]}"; do
        l=$(ls ${mapped}/${sample}*.out.bam)
        samtools merge --threads 4 ${merged}/${sample}.bam $l
    done
fi


# tag directories
if (($MERGING==1)); then
    for sample in "${samples[@]}"; do
        makeTagDirectory ${tagdir}/${sample} ${merged}/${sample}.bam -format bam -genome $genome
    done
else
    for sample in "${samples[@]}"; do
        makeTagDirectory ${tagdir}/${sample} ${mapped}/${sample}Aligned.sortedByCoord.out.bam -format bam -genome $genome
    done
fi


# analyze repeats
analyzeRepeats.pl rna hg19 -count genes -d IMR90-GroSeq > outputfile.txt


# DGE analysis
Rscript rbin/simple_DGE.R $INPUT_DIR $target_dir


# python script to output standard visuals for results

python3 script.py -i $resultsfile -o $OUTPUT_DIR

