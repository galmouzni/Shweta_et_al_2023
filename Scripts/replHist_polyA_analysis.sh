for xxx in {1..20}; do
    path_to_sample="/Volumes/Seb_GA/calcsub/slemaire/processed/D356T${xxx}.nsorted.bam"
    echo "BAM: ${path_to_sample}" >&2
    for ann in Data/ASF1_chaperone/annotations/variant_histone_and_histone_chaperone_genes.gff Data/ASF1_chaperone/annotations/replicative_histone_genes.gff Results/ASF1_chaperone/decile_expr_genes_increasing_syncCellsRNAseq/decile_9.gff; do
        path_to_annot=${ann}
        echo "Annot: ${path_to_annot}" >&2


        ## retrieve the reads
        alns="$( samtools view ${path_to_sample} | python3 ./Scripts/ASF1_chaperone/replHist_polyA_analysis_230929.py ${path_to_annot} )"
        if [ "$overwrite_sam" == "T" ]; then
            outbam="./Results/ASF1_chaperone/poly_adenylation_22T/$( basename ${ann} .gff )_D356T${xxx}_polyAreads.sam"
            echo "output BAM: ${outbam}" >&2
            echo "$alns" > ${outbam}
        fi

        res="$( echo "$alns" | wc -l )"
        echo D356T${xxx} $( basename ${ann} .gff ) ${res}
    done
done
