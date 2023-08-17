#!/bin/bash

## Usage: bash process.RNA.sh <run>
# <run>    accession number of the FASTQ file(s) in $PWD/data

#         i.e.  $PWD/data/<run>.fastq.gz (single-end)
#               $PWD/data/<run>.R[12].fastq.gz (paired-end)


set -e #exit script if any subcommand returns a non-zero status

echo "Args: ${@}"

source $PWD/scripts/env_bin.sh



run=$1


#set hisat2 arguments {-1 <m1> -2 <m2> | -U <r>}  
if [[ -f $PWD/data/$run.fastq.gz ]]; then
  
  #single-end reads
  args="-U $PWD/data/$run.fastq.gz"

elif [[ -f $PWD/data/$run.R1.fastq.gz && 
        -f $PWD/data/$run.R2.fastq.gz ]]; then
  
  #paired-end reads
  args="-1 $PWD/data/$run.R1.fastq.gz -2 $PWD/data/$run.R2.fastq.gz"

elif [[ -f $PWD/data/${run}_1.fastq.gz && 
        -f $PWD/data/${run}_2.fastq.gz ]]; then
  
  #paired-end reads
  args="-1 $PWD/data/${run}_1.fastq.gz -2 $PWD/data/${run}_2.fastq.gz"

else
    
  #exit if <run>.fastq.gz or <run>.R[12].fastq.gz not in $PWD/data/
  >&2 echo "$run(.R[12]).fastq.gz or $run(_[12]).fastq.gz not found in $PWD/data"; exit 1

fi

#output filename prefix
out=$PWD/processed/$run

#add suffix if specified
if [[ "x${@}" == x*"--suffix "* ]]
then
      out="${out}_$( echo "${@}" | awk -F '--suffix ' '{ print $2 }' | awk '{ print $1 }' )"
fi

#align and convert to BAM format
echo "$hisat2 -p 8 -x $idx $args | $samtools view -b - > $out.bam"

$hisat2 -p 8 -x $idx $args | $samtools view -b - > $out.bam

#sort BAM file (by read name)
echo "$samtools sort -n -m 2G -@ 8 $out.bam > $out.nsorted.bam"

$samtools sort -n -m 2G -@ 8 $out.bam > $out.nsorted.bam

#add tags required for markdup e.g. -m (mate score tag)
echo "$samtools fixmate -m -@ 8 $out.nsorted.bam $out.fixmate.bam"

$samtools fixmate -m -@ 8 $out.nsorted.bam $out.fixmate.bam

#sort BAM file (by position)
echo "$samtools sort -m 2G -@ 8 $out.fixmate.bam > $out.sorted.bam"

$samtools sort -m 2G -@ 8 $out.fixmate.bam > $out.sorted.bam
  
#mark duplicates
echo "$samtools markdup -@ 8 $out.sorted.bam $out.sorted.markdup.bam"    

$samtools markdup -@ 8 $out.sorted.bam $out.sorted.markdup.bam
   
#index BAM file
echo "$samtools index -b $out.sorted.markdup.bam"

$samtools index -b $out.sorted.markdup.bam

#remove temporary BAM files
echo "rm $out.bam $out.nsorted.bam $out.fixmate.bam $out.sorted.bam"

rm $out.bam $out.nsorted.bam $out.fixmate.bam $out.sorted.bam

#only primary alignments, excluding duplicates and reads with MAPQ lower than 2
opt="-Q 2 --primary --ignoreDup"

#paired-end option for read summarization
case "$args" in

  -1*) #count fragments
       opt="$opt -p" ;;
       
esac

#summarize gene-level counts (unstranded)
echo "$featureCounts -T 8 -t gene -O --minOverlap 5 -M -s 0 $opt -a $idx.gtf.gz -o ${out}_wholeGene.un $out.sorted.markdup.bam"
echo "$featureCounts -T 8 -t exon -O --minOverlap 5 -M -s 0 $opt -a $idx.gtf.gz -o ${out}_geneByExon.un $out.sorted.markdup.bam"
echo "$featureCounts -T 8 -t exon -f -O --minOverlap 5 -M -s 0 $opt -a $idx.gtf.gz -o ${out}_eachExon.un $out.sorted.markdup.bam"

$featureCounts -T 8 -t gene -O --minOverlap 5 -M -s 0 $opt -a $idx.gtf.gz -o ${out}_wholeGene.un $out.sorted.markdup.bam
$featureCounts -T 8 -t exon -O --minOverlap 5 -M -s 0 $opt -a $idx.gtf.gz -o ${out}_geneByExon.un $out.sorted.markdup.bam
$featureCounts -T 8 -t exon -f -O --minOverlap 5 -M -s 0 $opt -a $idx.gtf.gz -o ${out}_eachExon.un $out.sorted.markdup.bam

#summarize gene-level counts (stranded)
echo "$featureCounts -T 8 -t gene -O --minOverlap 5 -M -s 1 $opt -a $idx.gtf.gz -o ${out}_wholeGene.fw $out.sorted.markdup.bam"
echo "$featureCounts -T 8 -t exon -O --minOverlap 5 -M -s 1 $opt -a $idx.gtf.gz -o ${out}_geneByExon.fw $out.sorted.markdup.bam"
echo "$featureCounts -T 8 -t exon -f -O --minOverlap 5 -M -s 1 $opt -a $idx.gtf.gz -o ${out}_eachExon.fw $out.sorted.markdup.bam"

$featureCounts -T 8 -t gene -O --minOverlap 5 -M -s 1 $opt -a $idx.gtf.gz -o ${out}_wholeGene.fw $out.sorted.markdup.bam
$featureCounts -T 8 -t exon -O --minOverlap 5 -M -s 1 $opt -a $idx.gtf.gz -o ${out}_geneByExon.fw $out.sorted.markdup.bam
$featureCounts -T 8 -t exon -O --minOverlap 5 -M -f -s 1 $opt -a $idx.gtf.gz -o ${out}_eachExon.fw $out.sorted.markdup.bam

#summarize gene-level counts (reverse stranded)
echo "$featureCounts -T 8 -t gene -O --minOverlap 5 -M -s 2 $opt -a $idx.gtf.gz -o ${out}_wholeGene.rv $out.sorted.markdup.bam"
echo "$featureCounts -T 8 -t exon -O --minOverlap 5 -M -s 2 $opt -a $idx.gtf.gz -o ${out}_geneByExon.rv $out.sorted.markdup.bam"
echo "$featureCounts -T 8 -t exon -f -O --minOverlap 5 -M -s 2 $opt -a $idx.gtf.gz -o ${out}_eachExon.rv $out.sorted.markdup.bam"

$featureCounts -T 8 -t gene -O --minOverlap 5 -M -s 2 $opt -a $idx.gtf.gz -o ${out}_wholeGene.rv $out.sorted.markdup.bam
$featureCounts -T 8 -t exon -O --minOverlap 5 -M -s 2 $opt -a $idx.gtf.gz -o ${out}_geneByExon.rv $out.sorted.markdup.bam
$featureCounts -T 8 -t exon -f -O --minOverlap 5 -M -s 2 $opt -a $idx.gtf.gz -o ${out}_eachExon.rv $out.sorted.markdup.bam

#number of assigned reads (stranded)
n_fw=$(grep "Assigned" ${out}_wholeGene.fw.summary | cut -f 2)

#number of assigned reads (reverse stranded)
n_rv=$(grep "Assigned" ${out}_wholeGene.rv.summary | cut -f 2)

#keep strand-specific counts with higher number of assigned reads
echo "#assigned reads: $n_fw (forward) $n_rv (reverse)"

if (( $n_fw > $n_rv )); then echo "rm ${out}_*.rv*"; rm ${out}_*.rv*; fi

if (( $n_rv > $n_fw )); then echo "rm ${out}_*.fw*"; rm ${out}_*.fw*; fi


#message end of job
echo "JOB DONE" >&2

