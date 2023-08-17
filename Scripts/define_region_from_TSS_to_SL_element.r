#### Define region from TSS to SL element for the replicative histone genes

require(rtracklayer)
require(tidyverse)

rhbed <- rtracklayer::import.bed('Data/ASF1_chaperone/annotations/replicative_histone_genes_wHDE.sorted.bed')

rhbp <- './Data/ASF1_chaperone/annotations/replicative_histone_genes_wHDE.sorted.bed'
rhbed <- tibble(read.table(rhbp, h = F, sep = '\t'))
colnames(rhbed) <- c('ref', 'start_0.based', 'end', 'score', 'gene_name_38', 'strand')


slbp <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_stem_loop.bed'
slbed <- tibble(read.table(slbp, h = F, sep = '\t'))
colnames(slbed) <- c('ref_sl', 'start_0.based_sl', 'end_sl', 'score_sl', 'gene_name_38', 'strand_sl')

gvp <- './Data/GRCh38_GRCh37_gene_names_correspondences.tsv'
gvtab <- tibble(read.table(gvp, h = T, sep = '\t'))
colnames(gvtab) <- c('gene_id', 'gene_name_38', 'gene_name_37')


ttt <- left_join(
    rhbed,
    slbed,
    by = 'gene_name_38'
)
colnames(ttt)

ttt$tss_to_sl_start <- -1
ttt$tss_to_sl_end <- -1

poss <- ttt$strand == '+'
negs <- ttt$strand == '-'

ttt$tss_to_sl_start[poss] <- ttt$start_0.based[poss] + 1
ttt$tss_to_sl_end[negs] <- ttt$end[negs]

ttt$tss_to_sl_end[poss] <- ttt$start_0.based_sl[poss]
ttt$tss_to_sl_start[negs] <- ttt$end_sl[negs]

out <- ttt[c('ref', 'tss_to_sl_start', 'tss_to_sl_end', 'score', 'gene_name_38', 'strand')]
out <- left_join(
    out,
    gvtab,
    by = 'gene_name_38'
)
colnames(out) <- c('ref', 'start_1.based', 'end', 'score', 'gene_name_38', 'strand', 'gene_id', 'gene_name_37')
outp <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_TSS_to_stem_loop.tsv'
write.table(out, outp, sep = '\t', quote = F, col.names = F, row.names = F)

outbed <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_TSS_to_stem_loop.bed'
write.table(out[c('ref', 'start_1.based', 'end', 'score', 'gene_name_38', 'strand')], outbed, sep = '\t', quote = F, col.names = F, row.names = F)


####
