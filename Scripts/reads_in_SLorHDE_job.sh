
source $PWD/scripts/env_bin.sh

opt="-T 8 -O --minOverlap 1 -M -s 2 -Q 2 --primary --ignoreDup -p --countReadPairs"

fnb="${1}"
ann="${2}"

bam="$PWD/processed/D356T${fnb}.sorted.markdup.bam"
out="$PWD/processed/$( basename $ann .gtf )_D356T${fnb}.tsv"

$featureCounts -t gene -F GTF $opt -a $ann -o ${out} $bam

