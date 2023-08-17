
#### General Variables
bedtools=/Users/labo/bin/bedtools2/bin/bedtools
samtools=/Users/labo/bin/samtools-1.12/samtools

annots_hsl=Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_stem_loop.sorted.bed
annots_hde=Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_downstream_element.sorted.bed
annots_ups=Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_stem_loop_20bUpstream.sorted.bed

resdir=./Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/RNA-Seq_fragment_on_HDE_SL_and_upstream
mkdir -v $resdir
####

#### Some general parameters for commands
## for converting to fragments
#   # Set filter of properly aligned read pairs
filter="-f 2 -F 3840"

#   # position of paired-end fragments based on TLEN of upstream read (TLEN > 0) 
convert='( $3 ~ "^[0-9XY]" ) && ( $9 > 0 ) { pos = $4 - 1
                                             end = pos + $9     
                                             print $3, pos, end }'
####


#### MAIN
for xxx in $( seq 5 20 ); do
    inbam=/Volumes/Seb_GA/calcsub/slemaire/processed/D356T${xxx}.sorted.markdup.bam
    ls $inbam

    #   # final destination
    out=${resdir}/$( basename $inbam .bam )

    #   # filter for reads with minimum mapping quality and without jumps in mapping
    $samtools view -h -q 10 $inbam | awk '( $1 ~ "^@" || $6 !~ "N" ) { print $0 }' OFS='\t' | samtools view -b > $out.mapq10.Nfiltered.bam && $samtools index $out.mapq10.Nfiltered.bam #&

    #   # convert into fragment annotations, with a maximum length, and filter for fragments overlapping the replicative histone stem loop, the position 20b upstream of the stem loop and the histone downstream element, the three together
    $samtools view $out.mapq10.Nfiltered.bam | awk "$convert" OFS="\t" | awk '{ if ( $3 - $2 <= 500 ) { print $0 } }' OFS="\t" | $bedtools intersect -u -a stdin -b ${annots_hsl} | $bedtools intersect -u -a stdin -b ${annots_ups} | $bedtools intersect -u -a stdin -b ${annots_hde} | gzip > $out.mapq10.Nfiltered.hSL.20upstrOfSL.HDE.frag500bMax.bed.gz #&

    #   # counts the HDE-overlapping fragments per replicative histone gene (one HDE per gene)
    $bedtools intersect -c -a ${annots_hde} -b $out.mapq10.Nfiltered.hSL.20upstrOfSL.HDE.frag500bMax.bed.gz | gzip > $out.mapq10.Nfiltered.hSL.20upstrOfSL.HDE.frag500bMax.countPerHDE.bed.gz
done



####
