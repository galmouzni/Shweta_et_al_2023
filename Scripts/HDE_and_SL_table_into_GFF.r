
if (TRUE) {
    fpv <- c("Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_stem_loop.bed", "Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_downstream_element.bed")
    for (fpath in fpv) {
        tab <- read.table(
            fpath,
            h=F,
            stringsAsFactors = F,
            sep = '\t'
        )
        colnames(tab) <- c('chrom', 'start_0.based', 'end', 'score', 'histone', 'strand')

        gr <- GRanges(
            seqnames = Rle(as.character(tab$chrom)),
            ranges = IRanges(start = tab$start_0.based + 1, end = tab$end),
            strand = tab$strand,
            type = 'gene', #ifelse(grepl('stem_loop', fpath), 'stem_loop', 'histone_downstream_element'),
            gene_id = tab$histone
        )

        opath <- sub('\\.bed$', '.gff', fpath)
        export.gff(con = opath, object = gr, format = "gff3")
    }
}
