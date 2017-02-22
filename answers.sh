#! /usr/bin/env bash 

datasets='/Users/lorberbd/data-sets/'


# Use BEDtools intersect to identify the size of the largest overlap
# between CTCF and H3K4me3 locations. 

h3k4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"
TFBS="$datasets/bed/encode.tfbs.chr22.bed.gz"

answer_1=$(bedtools intersect -a $h3k4me3 -b $TFBS \
    | awk 'BEGIN {OFS="\t"} {print  $0, $3- $2}'  \
    | sort -k11nr \
    | head -n1 \
    | cut -f11)

echo "answer-1: $answer_1"

#Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#19,000,500 on chr22 of hg19 genome build. Report the GC content as a
#fraction (e.g., 0.50). 

region="$datasets/bed/q2.hg19.bed.gz"
hg19="$datasets/fasta/hg19.chr22.fa"

answer_2=$(bedtools nuc -fi $hg19 -bed $region \
    | cut -f5 \
    | tail -n1)

echo "answer-2: $answer_2"


#Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

#First I generated a bed file that contains CTCF binding sites from the
#ENCODE TFBS file we had and called it $ctcfpeaks below. Here is the code
# I used: $ gzcat encode.tfbs.chr22.bed.gz | grep -w "CTCF" >
# ctcfpeaks.chr22.bed.gz

#I'm sure there is an easier way to do this, but
#I think this works too..maybe? 

ctcfpeaks="$datasets/bed/ctcfpeaks.chr22.bed.gz"
ctcf="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_3=$(bedtools map -c 4 -o mean -a $ctcfpeaks -b $ctcf \
    | sort -k5nr \
    | head -n1 \
    | cut -f1-3)

echo "answer-3: $answer_3"

#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
#of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz. Report
#the gene name (e.g., 'ABC123')


tss="$datasets/bed/tss.hg19.chr22.bed.gz"
hg19="$datasets/genome/hg19.genome"
ctcfhela="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_4=$(bedtools slop -l 1000 -r 0 -s -i $tss -g $hg19 \
    | sort \
    | bedtools map -c 4 -o median -a - -b $ctcfhela \
    | sort -k7nr \
    | head -n1 \
    | cut -f4)
echo "answer-4: $answer_4"

#Use BEDtools to identify the longest interval on chr22 that is not
#covered by genes.hg19.bed.gz. Report the interval like chr1:100-500

genes="$datasets/bed/genes.hg19.bed.gz"
hg19genome="$datasets/genome/hg19.genome"

answer_5=$(bedtools complement -i $genes -g $hg19genome \
    | awk 'BEGIN {OFS="\t"} ($1 == "chr22") {print $1, $2, $3, $3-$2}' \
    | sort -k4nr \
    | cut -f1-3 \
    | head -n1)
echo "answer-5: $answer_5"

