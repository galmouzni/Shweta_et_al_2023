
source $PWD/scripts/env_bin.sh

opt="-T 8 -O --minOverlap 1 -M -s 2 -Q 2 --primary --ignoreDup -p --countReadPairs"

spl="${1}" # ex: D356T1, D1182-1176T03, ... [name or id of the sample BAM file (<spl>.bam)]
ann="${2}" # ex: ${path}/annot.gtf  [path to annotations table, GTF format]
anntype="${3}" # ex: gene, CDS, exon, ...  [annotation types in the annotations table to consider]

bam="$PWD/processed/${spl}.sorted.markdup.bam"
out="$PWD/processed/$( basename $ann .gtf )_${spl}.tsv"

$featureCounts -t $anntype -F GTF $opt -a $ann -o ${out} $bam

