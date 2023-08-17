for xxx in {1..20}; do
    path_to_sample="/Volumes/Seb_GA/calcsub/slemaire/processed/D356T${xxx}.nsorted.bam"
    for ann in Data/ASF1_chaperone/annotations/variant_histone_and_histone_chaperone_genes.gff Data/ASF1_chaperone/annotations/replicative_histone_genes.gff Results/ASF1_chaperone/decile_expr_genes_increasing_syncCellsRNAseq/decile_9.gff; do
        path_to_annot=${ann}
        res="$( samtools view ${path_to_sample} | python3 ./Scripts/ASF1_chaperone/replHist_polyA_analysis_230630.py ${path_to_annot} | wc -l )"
        echo D356T${xxx} $( basename ${ann} .gff ) ${res}
    done
done