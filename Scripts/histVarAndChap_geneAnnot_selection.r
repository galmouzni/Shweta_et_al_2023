
require(rtracklayer)
require(tidyverse)

genes <- rtracklayer::import.gff('./Data/GRCh38.104.sorted.geneOnly.gtf')

hvct <- read.table('./Data/GtexTcgaGeneExpr/histone_chaperone_complete_list.csv', h=T, sep='\t')

hvcg <- genes[(genes$gene_id %in% hvct$EnsemblGeneId)]

is_replicative_v <- left_join(
    as.data.frame(mcols(hvcg)),
    hvct[c('EnsemblGeneId', 'subClass')] %>% mutate('gene_id'=EnsemblGeneId),
    by = 'gene_id'
)$subClass == 'replicative'

hvcg <- hvcg[!is_replicative_v]

rtracklayer::export.gff3(hvcg, './Data/ASF1_chaperone/annotations/variant_histone_and_histone_chaperone_genes.gff', )

####
