
{
    require(tidyverse)
    require(rtracklayer)
    require(plotly)
}

## General data
if (TRUE) {
    data_dir <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/RNA-Seq_fragment_on_HDE_SL_and_upstream'
    read_count_dir <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/featureCounts/on_coding_sequence_of_replicative_histone_genes' #'/Volumes/Seb_GA/calcsub/slemaire/processed'
    figdir <- paste0("./Images/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE"); dir.create(figdir, showWarnings = F)

    all_annot <- rtracklayer::import.gff('./Data/GRCh38.104.sorted.gtf')
    genes <- all_annot[all_annot$type == "gene"]

    ## get EnsEMBL annotation of coding sequenes for replicative histones genes
    cds_rhc_p <- './Data/ASF1_chaperone/annotations/replicative_histone_CDS.gtf'
    cdsgr <- rtracklayer::import.gff(cds_rhc_p)

    #### Extract promoter regions of the replicative histone genes
    if (FALSE) {
        resdir <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE'
        
        hc.genes <- genes[genes$gene_name %in% cdsgr$gene_name]

        hc.p1k <- hc.genes
        mcols(hc.p1k) <- mcols(hc.p1k)[c('gene_id', 'gene_name')]
        log_pos_strand <- strand(hc.p1k) == '+'
        log_neg_strand <- !log_pos_strand
        
        end(hc.p1k[log_pos_strand]) <- start(hc.p1k[log_pos_strand]) - 1
        start(hc.p1k[log_pos_strand]) <- end(hc.p1k[log_pos_strand]) - 1000 + 1
        start(hc.p1k[log_neg_strand]) <- end(hc.p1k[log_neg_strand]) + 1
        end(hc.p1k[log_neg_strand]) <- start(hc.p1k[log_neg_strand]) + 1000 - 1

        rtracklayer::export.gff(hc.p1k, file.path(resdir, 'histone_promoter1kbp.gtf'))

        hc.p1kg2h <- hc.p1k
        end(hc.p1kg2h[log_pos_strand]) <- end(hc.p1kg2h[log_pos_strand]) + 200
        start(hc.p1kg2h[log_neg_strand]) <- start(hc.p1kg2h[log_neg_strand]) - 200

        rtracklayer::export.gff(hc.p1kg2h, file.path(resdir, 'histone_promoter1kbpGene200bp.gtf'))

    }
    ####

    spl_idx_vec <- c(
        "1" = "siCtrl_0h_1",
        "2" = "siCtrl_0h_2",
        "3" = "siCtrl_2h_1",
        "4" = "siCtrl_2h_2",
        "5" = "siCtrl_5h_1",
        "6" = "siCtrl_5h_2",
        "7" = "siAsf1_0h_1",
        "8" = "siAsf1_0h_2",
        "9" = "siAsf1_2h_1",
        "10" = "siAsf1_2h_2",
        "11" = "siAsf1_5h_1",
        "12" = "siAsf1_5h_2",
        "13" = "siCtrl_asyn_1",
        "14" = "siCtrl_asyn_2",
        "15" = "siAsf1_asyn_1",
        "16" = "siAsf1_asyn_2",
        "17" = "siCtrl_nasc_1",
        "18" = "siAsf1_nasc_1",
        "19" = "siCtrl_nasc_2",
        "20" = "siAsf1_nasc_2"
    )
    spl_idx_vec_rev <- as.integer(names(spl_idx_vec))
    names(spl_idx_vec_rev) <- as.vector(spl_idx_vec)
    spl_order <- as.vector(spl_idx_vec)[c(1:16, 17, 19, 18, 20)]


    spl_replicates_list <- list(
        "siCtrl_0h" = c(1, 2), # c("siCtrl_0h_1", "siCtrl_0h_2"),
        "siCtrl_2h" = c(3, 4), # c("siCtrl_2h_1", "siCtrl_2h_2"),
        "siCtrl_5h" = c(5, 6), # c("siCtrl_5h_1", "siCtrl_5h_2"),
        "siAsf1_0h" = c(7, 8), # c("siAsf1_0h_1", "siAsf1_0h_2"),
        "siAsf1_2h" = c(9, 10), # c("siAsf1_2h_1", "siAsf1_2h_2"),
        "siAsf1_5h" = c(11, 12), # c("siAsf1_5h_1", "siAsf1_5h_2"),
        "siCtrl_asyn" = c(13, 14), # c("siCtrl_asyn_1", "siCtrl_asyn_2"),
        "siAsf1_asyn" = c(15, 16), # c("siAsf1_asyn_1", "siAsf1_asyn_2")
        "siCtrl_nasc" = c(17, 19), # c("siCtrl_nasc_1", "siCtrl_nasc_2"),
        "siAsf1_nasc" = c(18, 20) # c("siAsf1_nasc_1", "siAsf1_nasc_2")
    )


    ctrl_to_sirna_relation <- list(
        "0h" = c("siCtrl_0h", "siAsf1_0h"),
        "2h" = c("siCtrl_2h", "siAsf1_2h"),
        "5h" = c("siCtrl_5h", "siAsf1_5h"),
        "asyn" = c("siCtrl_asyn", "siAsf1_asyn"),
        "nasc" = c("siCtrl_nasc", "siAsf1_nasc")
    )


    desdf <- tibble(
        'sample' = factor(paste0('D356T', 1:20)),
        'treatment' = factor(c("siCtrl", "siCtrl", "siCtrl", "siCtrl", "siCtrl", "siCtrl", "siAsf1", "siAsf1", "siAsf1", "siAsf1", "siAsf1", "siAsf1", "siCtrl", "siCtrl", "siAsf1", "siAsf1", "siCtrl", "siAsf1", "siCtrl", "siAsf1"), levels = c('siCtrl', 'siAsf1'), ordered = TRUE),
        'time_rel' = factor(c(0, 0, 2, 2, 5, 5, 0, 0, 2, 2, 5, 5, 'asyn', 'asyn', 'asyn', 'asyn', 'nasc', 'nasc', 'nasc', 'nasc'), levels = c(0, 2, 5, 'asyn', 'nasc')),
        'replicate' = factor(c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 2))
    )

    treat_colpal <- c('siCtrl' = '#3174a1', 'siAsf1' = '#c03d3e')
    sufSplTypel <- list('asyn' = 'asynchroCells', 'synchro' = 'acrossSphase', 'nasc' = 'inNascentRNA')
}


## Compile the counts of unclived transcripts (~fragments)
ctab <- list()
for (filename in list.files(data_dir)) {
    if (grepl('D356T[0-9]*\\.sorted\\.markdup\\.mapq10\\.Nfiltered\\.hSL\\.20upstrOfSL\\.HDE\\.frag500bMax\\.countPerHDE\\.bed\\.gz', filename)) {
        message(filename)
        intab <- read.table(paste0(data_dir, '/', filename))
        colnames(intab) <- c('ref', 'pos_0.based', 'end', 'trash', 'repl_histone', 'strand', 'count')
        ctab[[filename]] <- intab
    }
}

# lapply(ctab, function(xxx) {
#     head(xxx)
# })

gdata <- tibble(ctab[[1]])
gdata <- gdata[1:ncol(gdata) - 1]

for (name in names(ctab)) {
    gdata[[strsplit(name, '\\.')[[1]][1]]] <- ctab[[name]][[7]]
}

## normalize the counts by clivage window length
gdata_lenNorm <- gdata
for (coln in colnames(gdata_lenNorm)[grepl('^D356T', colnames(gdata_lenNorm))]) {
    gdata_lenNorm[[coln]] <- gdata_lenNorm[[coln]] / 49 #gdata_lenNorm$Length
}


## Compile the counts of reads per replicative histone
rtab <- list()
seq_depth_list <- list()
for (filename in list.files(read_count_dir)) {
    # if (grepl('D356T[0-9]*_wholeGene.rv$', filename)) {
    if (grepl('replicative_histone_CDS_D356T[0-9]*.tsv$', filename)) {
        message(filename)
        intab <- read.table(paste0(read_count_dir, '/', filename), h = T)
        seq_depth_list[[filename]] <- sum(intab[[7]])
        intab <- intab[intab$Geneid %in% cdsgr$gene_id,]
        rtab[[filename]] <- intab
    }
}

# dim(rtab[[1]])

#### Compile the total mRNA read counts for the replication histone genes
if (TRUE) {
    ## reduce the count tables in one unique table
    rcdata <- tibble(rtab[[1]])
    rcdata <- rcdata[1:ncol(rcdata) - 1]
    rcdata['repl_histone'] <- left_join(
        rcdata['Geneid'] %>% mutate(gene_id = Geneid),
        as.data.frame(mcols(genes))[c('gene_id', 'gene_name')],
        by = 'gene_id'
    )[['gene_name']]


    for (name in names(rtab)) {
        spl_name <- strsplit(strsplit(name, '\\.')[[1]][1], '_')[[1]][4]
        rcdata[[spl_name]] <- rtab[[name]][[7]]
    }

    ## give same order of replicative histone than for unprocessed transcript table
    rcdata <- left_join(
        gdata['repl_histone'],
        rcdata,
        by = 'repl_histone'
    )
    rcdata <- rcdata[colnames(rcdata)[c(2:7, 1, 8:ncol(rcdata))]]

    ## look for the distribution of the repl. histone for their gene expression
    if (TRUE) {
        figdata <- rcdata[grepl('^(repl_histone|D356T[0-9]*)$', colnames(rcdata))]
        
        # ## normalize by total sequencing depth
        # for (name in names(rtab)) {
        #     spl_name <- strsplit(strsplit(name, '\\.')[[1]][1], '_')[[1]][1]
        #     figdata[[spl_name]] <- figdata[[spl_name]] / ( seq_depth_list[[name]] / 1000 )
        # }

        figdata <- pivot_longer(figdata, colnames(figdata)[2:ncol(figdata)], names_to = "condition", values_to = "expression")

        figdata <- figdata[order(figdata$condition),]
        figdata[['rank_by_condition']] <- unlist(tapply(figdata$expression, figdata$condition, function(xxx) {order(order(xxx))}))

        figdata <- figdata[order(figdata$repl_histone),]
        figdata[['rank_by_repl_histone']] <- unlist(tapply(figdata$expression, figdata$repl_histone, function(xxx) {order(order(xxx))}))

        ## prevent null values
        figdata$expression[figdata$expression == 0] <- 0.5

        ## build expression repartition by condition
        if (TRUE) {
            figdata <- figdata[order(figdata$rank_by_condition),]
            fig <- ggplot(
                figdata,
                aes(
                    x = expression,
                    y = rank_by_condition,
                    col = condition
                )
            ) +
                geom_line() +
                scale_x_continuous(trans = 'log10')
            
            figname <- paste0(figdir, '/', 'repl_histone_expression_repartition_by_condition.pdf'); message(figname)
            pdf(figname)
            print(fig)
            bouh <- dev.off()
        }

        ## build expression repartition by replicative histone
        if (TRUE) {
            figdata <- figdata[order(figdata$rank_by_repl_histone),]
            fig <- ggplot(
                figdata,
                aes(
                    x = expression,
                    y = rank_by_repl_histone,
                    col = repl_histone
                )
            ) +
                geom_line() +
                geom_vline(xintercept=45) +
                scale_x_continuous(trans = 'log10')
            
            figname <- paste0(figdir, '/', 'repl_histone_expression_repartition_by_replHistone.pdf'); message(figname)
            pdf(figname)
            print(fig)
            bouh <- dev.off()
        }

    }

    ## normalize the counts by annotation length
    rcdata_lenNorm <- rcdata
    for (coln in colnames(rcdata_lenNorm)[grepl('^D356T', colnames(rcdata_lenNorm))]) {
        rcdata_lenNorm[[coln]] <- rcdata_lenNorm[[coln]] / rcdata_lenNorm$Length
    }
}

## determine replicative histone genes with enough reads
min_expr <- 45 #100
min_expr_log <- apply(rcdata[8:ncol(rcdata)], 1, function(xxx) {
    all(xxx > min_expr)
})


#### Normalizes expression and boxplot separating total and unprocessed for total RNA-Seq experiment
if (TRUE) {
    ## Plot expression through S-phase without normalization between the conditions
    if (TRUE) {
        figdata <- rcdata[min_expr_log,]
        figdata <- pivot_longer(figdata[7:ncol(figdata)], colnames(figdata)[8:ncol(figdata)], names_to = 'sample', values_to = 'expression')
        # figdata <- figdata[!grepl('T(17|18|19|20)$', figdata$sample),]
        figdata <- left_join(
            figdata,
            desdf,
            by = 'sample'
        )
        figdata$group <- paste0(figdata$treatment, '_', figdata$time_rel, '_', figdata$replicate)
        figdata$group <- factor(figdata$group, levels = c('siCtrl_0_1', 'siCtrl_0_2', 'siCtrl_2_1', 'siCtrl_2_2', 'siCtrl_5_1', 'siCtrl_5_2', 'siCtrl_asyn_1', 'siCtrl_asyn_2', 'siCtrl_nasc_1', 'siCtrl_nasc_2', 'siAsf1_0_1', 'siAsf1_0_2', 'siAsf1_2_1', 'siAsf1_2_2', 'siAsf1_5_1', 'siAsf1_5_2', 'siAsf1_asyn_1', 'siAsf1_asyn_2', 'siAsf1_nasc_1', 'siAsf1_nasc_2'), ordered = TRUE)

        fig <- ggplot(
            figdata,
            aes(
                x = time_rel,
                y = expression,
                group = group,
                fill = treatment
            )
        ) +
            geom_boxplot() +
            scale_y_continuous(trans = 'log10') +
            scale_fill_manual(values = treat_colpal)
        
        figname <- paste0(figdir, '/', 'Total_replicationHistone_raw_expression_acrossSphase_bxplt', '.pdf'); message(figname)
        pdf(figname)
        print(fig)
        bouh <- dev.off()
    }

    ## load DESeq2 results (for getting normalization factors)
    if (TRUE) {
        cat.desq.path <- './Results/ASF1_chaperone/DESeq2_results/CtrlVsSiASF1_while_AsyncVsSphase_effect_DESeq2Results.RDS'
        cat.dds <- readRDS(cat.desq.path)
        cat.dds.ld <- tibble(as.data.frame(cat.dds@colData@listData))

        nas.desq.path <- './Results/ASF1_chaperone/DESeq2_results/CtrlVsSiASF1_inNascentRNAs_DESeq2Results.RDS'
        nas.dds <- readRDS(nas.desq.path)
        nas.dds.ld <- tibble(as.data.frame(nas.dds@colData@listData))
    }
    
    ## normalize expression by DESeq2 normalization factors
    if (TRUE) {
        ngdata_lenNorm <- gdata_lenNorm
        nrcdata_lenNorm <- rcdata_lenNorm

        for (coln in cat.dds.ld$sample) {
            ngdata_lenNorm[[coln]] <- ngdata_lenNorm[[coln]] / cat.dds.ld$sizeFactor[cat.dds.ld$sample == coln]
            nrcdata_lenNorm[[coln]] <- nrcdata_lenNorm[[coln]] / cat.dds.ld$sizeFactor[cat.dds.ld$sample == coln]
        }

        for (coln in nas.dds.ld$sample) {
            ngdata_lenNorm[[coln]] <- ngdata_lenNorm[[coln]] / nas.dds.ld$sizeFactor[nas.dds.ld$sample == coln]
            nrcdata_lenNorm[[coln]] <- nrcdata_lenNorm[[coln]] / nas.dds.ld$sizeFactor[nas.dds.ld$sample == coln]
        }

        ngdata_lenNorm <- ngdata_lenNorm[c(5, 7:ncol(ngdata_lenNorm))]
        nrcdata_lenNorm <- nrcdata_lenNorm[7:ncol(nrcdata_lenNorm)]
    }

    ## Plot normalized expression of replicative histones through S-phase
    if (TRUE) {
        figdata <- nrcdata_lenNorm[min_expr_log,]
        figdata <- pivot_longer(figdata, colnames(figdata)[2:ncol(figdata)], names_to = 'sample', values_to = 'expression')
        # figdata <- figdata[!grepl('T(17|18|19|20)$', figdata$sample),]
        figdata <- left_join(
            figdata,
            desdf,
            by = 'sample'
        )
        figdata$group <- paste0(figdata$treatment, '_', figdata$time_rel, '_', figdata$replicate)
        figdata$group <- factor(figdata$group, levels = c('siCtrl_0_1', 'siCtrl_0_2', 'siCtrl_2_1', 'siCtrl_2_2', 'siCtrl_5_1', 'siCtrl_5_2', 'siCtrl_asyn_1', 'siCtrl_asyn_2', 'siCtrl_nasc_1', 'siCtrl_nasc_2', 'siAsf1_0_1', 'siAsf1_0_2', 'siAsf1_2_1', 'siAsf1_2_2', 'siAsf1_5_1', 'siAsf1_5_2', 'siAsf1_asyn_1', 'siAsf1_asyn_2', 'siAsf1_nasc_1', 'siAsf1_nasc_2'), ordered = TRUE)

        fig <- ggplot(
            figdata,
            aes(
                x = time_rel,
                y = expression,
                group = group,
                fill = treatment
            )
        ) +
            geom_boxplot() +
            scale_y_continuous(trans = 'log10') +
            scale_fill_manual(values = treat_colpal)
        
        figname <- paste0(figdir, '/', 'Total_replicationHistone_normalized_expression_acrossSphase_bxplt', '.pdf'); message(figname)
        pdf(figname)
        print(fig)
        bouh <- dev.off()
    }
}


#### Look at the relation between total expression and amount of unprocessed transcripts (2D total vs unprocessed)
if (TRUE) {
    figdata <- bind_rows(
        ngdata_lenNorm[min_expr_log,] %>% pivot_longer(colnames(ngdata_lenNorm)[2:ncol(ngdata_lenNorm)], names_to = 'sample', values_to = 'expression') %>% mutate(val_type = 'unprocessed'),
        nrcdata_lenNorm[min_expr_log,] %>% pivot_longer(colnames(ngdata_lenNorm)[2:ncol(ngdata_lenNorm)], names_to = 'sample', values_to = 'expression') %>% mutate(val_type = 'whole_gene')
    )

    figdata <- figdata[!grepl('T(17|18|19|20)$', figdata$sample),]
    figdata$expression[figdata$expression == 0] <- 0.5
    figdata <- pivot_wider(figdata, c('repl_histone', 'sample'), names_from = val_type, values_from = expression)
    # figdata[(is.na(figdata$whole_gene)),]

    figdata <- left_join(
        figdata,
        desdf,
        by = 'sample'
    )

    figdata <- (figdata %>% mutate(synchro = c('synchro', 'asyn', 'nasc')[(1 + (time_rel %in% c('asyn', 'nasc')) + (time_rel == 'nasc'))]))

    # 2D plot Unprocessed ~ f(Total)
    for (synchro in unique(figdata$synchro)) {
        # synchro <- "synchro"

        fig <- ggplot(
            figdata[figdata$synchro == synchro,],
            aes(
                x = log10(whole_gene),
                y = log10(unprocessed),
                # shape = treatment,
            )
        )
        # if (synchro == 'synchro') {
        #     fig <- fig +
        #         facet_wrap(~treatment, nrow = 1) +
        #         geom_point(aes(col = time_rel))
        # } else {
        #     fig <- fig +
        #         geom_point(aes(col = treatment)) +
        #         scale_color_manual(values = treat_colpal)
        # }
        fig <- fig +
            facet_wrap(~time_rel, nrow = 1) +
            geom_point(aes(col = treatment)) +
            scale_color_manual(values = treat_colpal)
        fig <- fig +
            # scale_x_continuous(trans = 'log10') +
            # scale_y_continuous(trans = 'log10') +
            xlim(log10(min(figdata$whole_gene)), log10(max(figdata$whole_gene))) +
            ylim(log10(min(figdata$unprocessed)), log10(max(figdata$unprocessed))) +
            theme(text = element_text(size = 20))
        
        figname <- paste0(figdir, '/', 'TotalvsUnprocessed_replicationHistone_normalized_expression_', sufSplTypel[[synchro]], '_dot', '.pdf'); message(figname)
        pdf(figname, width = 7 * c(1.5, 2.5)[1+(synchro == 'synchro')])
        print(fig)
        bouh <- dev.off()
    }

    # boxplot about Unprocessed quantities
    for (synchro in unique(figdata$synchro)) {
        fig <- ggplot(
            figdata[figdata$synchro == synchro,],
            aes(
                x = time_rel,
                y = log10(unprocessed),
                linetype = replicate,
                group = factor(paste0(treatment, '_', time_rel, '_', replicate), levels = c('siCtrl_0_1', 'siCtrl_0_2', 'siCtrl_2_1', 'siCtrl_2_2', 'siCtrl_5_1', 'siCtrl_5_2', 'siCtrl_asyn_1', 'siCtrl_asyn_2', 'siCtrl_nasc_1', 'siCtrl_nasc_2', 'siAsf1_0_1', 'siAsf1_0_2', 'siAsf1_2_1', 'siAsf1_2_2', 'siAsf1_5_1', 'siAsf1_5_2', 'siAsf1_asyn_1', 'siAsf1_asyn_2', 'siAsf1_nasc_1', 'siAsf1_nasc_2'), ordered = TRUE)
            )
        )
        fig <- fig +
            geom_boxplot(aes(fill = treatment)) +
            scale_fill_manual(values = treat_colpal)
        fig <- fig +
            # scale_x_continuous(trans = 'log10') +
            # scale_y_continuous(trans = 'log10') +
            ylim(log10(min(figdata$unprocessed)), log10(max(figdata$unprocessed))) +
            theme(text = element_text(size = 20))
        
        figname <- paste0(figdir, '/', 'Unprocessed_replicationHistone_normalized_expression_', sufSplTypel[[synchro]], '_bxplt', '.pdf'); message(figname)
        pdf(figname, width = 7 * c(0.7, 1)[1+(synchro == 'synchro')])
        print(fig)
        bouh <- dev.off()
    }

    # boxplot about Total expression
    for (synchro in unique(figdata$synchro)) {
        fig <- ggplot(
            figdata[figdata$synchro == synchro,],
            aes(
                x = time_rel,
                y = log10(whole_gene),
                linetype = replicate,
                group = factor(paste0(treatment, '_', time_rel, '_', replicate), levels = c('siCtrl_0_1', 'siCtrl_0_2', 'siCtrl_2_1', 'siCtrl_2_2', 'siCtrl_5_1', 'siCtrl_5_2', 'siCtrl_asyn_1', 'siCtrl_asyn_2', 'siCtrl_nasc_1', 'siCtrl_nasc_2', 'siAsf1_0_1', 'siAsf1_0_2', 'siAsf1_2_1', 'siAsf1_2_2', 'siAsf1_5_1', 'siAsf1_5_2', 'siAsf1_asyn_1', 'siAsf1_asyn_2', 'siAsf1_nasc_1', 'siAsf1_nasc_2'), ordered = TRUE)
            )
        )
        fig <- fig +
            geom_boxplot(aes(fill = treatment)) +
            scale_fill_manual(values = treat_colpal)
        fig <- fig +
            # scale_x_continuous(trans = 'log10') +
            # scale_y_continuous(trans = 'log10') +
            ylim(log10(min(figdata$whole_gene)), log10(max(figdata$whole_gene))) +
            theme(text = element_text(size = 20))
        
        figname <- paste0(figdir, '/', 'WholeGene_replicationHistone_normalized_expression_', sufSplTypel[[synchro]], '_bxplt', '.pdf'); message(figname)
        pdf(figname, width = 7 * c(0.7, 1)[1+(synchro == 'synchro')])
        print(fig)
        bouh <- dev.off()
    }
}


#### Compute the fractions of unclived mRNA (in all condition and experiment type [total or nascent RNA])
if (TRUE) {
    ## check tables of total expression and of unprocessed transcript have the same order of replicative histone genes
    if (!all(gdata_lenNorm$repl_histone == rcdata_lenNorm$repl_histone)) { 
        message('!!! Different order between "gdata" and "rcdata". Cannot go further !!!')
    } else {
        ## compute fraction by sample
        fdata <- gdata_lenNorm
        for (nameid in names(spl_idx_vec)) {
            spl_name <- paste0("D356T", nameid)
            # spl_name <- strsplit(name, '\\.')[[1]][1]
            fdata[[spl_name]] <- gdata_lenNorm[[spl_name]] / rcdata_lenNorm[[spl_name]]
        }

        ## average on the replicates
        afdata <- fdata[1:6]
        for (condname in names(spl_idx_vec_rev)) {
            afdata[[condname]] <- fdata[[paste0("D356T", spl_idx_vec_rev[condname])]]
        }

        mafdata <- fdata[1:6]
        for (condname in names(spl_replicates_list)) {
            mafdata[[condname]] <- apply(fdata[paste0("D356T", spl_replicates_list[[condname]])], 1, mean, na.rm = T)
        }

        ## establish stable set of replication histone (with valid fraction in all conditions)
        if (TRUE) {
            rhlog <- rep(TRUE, length(afdata$repl_histone))
            for (syntime in names(ctrl_to_sirna_relation)) {
                for (repli in 1:2) {
                    xcol <- paste0("siCtrl_", syntime, "_", repli)
                    ycol <- paste0("siAsf1_", syntime, "_", repli)

                    clog <- apply(afdata[c(xcol, ycol)], 1, function(xxx) {
                        !any(is.nan(xxx) | is.infinite(xxx))
                    })
                    # print(clog)
                    rhlog <- rhlog & clog
                }
                # xcol <- paste0("siCtrl_", syntime)
                # ycol <- paste0("siAsf1_", syntime)

                # clog <- apply(afdata[c(xcol, ycol)], 1, function(xxx) {
                #     !any(is.nan(xxx) | is.infinite(xxx))
                # })
                # print(clog)
                rhlog <- rhlog & clog
            }
        }

        row_log <- rhlog & min_expr_log & afdata$repl_histone != 'H1-2'
        # sum(row_log)

        ## plot total expression vs unprocessed across time
        if (TRUE) {
            ## the plot itself
            if (TRUE) {
                figdata <- afdata[row_log, c('repl_histone', colnames(afdata)[grepl('^si(Ctrl|Asf1)_', colnames(afdata))])]
                figdata <- pivot_longer(figdata, colnames(figdata)[2:ncol(figdata)], names_to = 'condition', values_to = 'rel_unpr')
                figdata <- bind_cols(
                    figdata,
                    (tibble(data.frame(do.call(rbind, strsplit(figdata$condition, split='_')))) %>% mutate(treatment = X1, rel_time = X2, replicate = X3))[c('treatment', 'rel_time', 'replicate')]
                )
                figdata$treatment <- factor(figdata$treatment, levels = c('siCtrl', 'siAsf1'), ordered = TRUE)
                # all(figdata$condition == paste0(figdata$treatment, '_', figdata$rel_time))
                figdata$condition <- factor(figdata$condition, levels = spl_order)
                figdata$facetx <- c('synchro', 'asyn', 'nasc')[(1 + (figdata$rel_time %in% c('asyn', 'nasc')) + (figdata$rel_time == 'nasc'))]

                for (synchro in unique(figdata$facetx)) {
                    fig <- ggplot(
                        figdata[figdata$facetx == synchro,],
                        aes(
                            x = rel_time,
                            y = rel_unpr,
                            linetype = replicate,
                            group = condition,
                            fill = treatment
                        )
                    ) +
                        geom_boxplot() +
                        scale_fill_manual(values = treat_colpal) +
                        ylim(0, ifelse(synchro == 'nasc', max(figdata$rel_unpr[figdata$rel_time == "nasc"]), max(figdata$rel_unpr[figdata$rel_time != "nasc"]))) +
                        theme(text = element_text(size = 20))
                    
                    figname <- paste0(figdir, '/', 'relAmountUnprocessed_through_time_', sufSplTypel[[synchro]], '_boxplot.pdf'); message(figname)
                    pdf(figname, width = 7 * c(0.7, 1)[1 + (synchro == "synchro")])
                    print(fig)
                    bouh <- dev.off()

                    ## associated statistics
                    {
                        if (synchro %in% c('nasc', 'asyn')) {
                            message(paste0('Treatment effect on ', synchro, ' cells.'))

                            glmres <- glm(
                                data = figdata[figdata$facetx == synchro,],
                                formula = rel_unpr ~ replicate + treatment + treatment:replicate
                            )
                            print(summary(glmres))

                            # t.test(
                            #     x = figdata$rel_unpr[figdata$facetx == synchro & figdata$treatment == 'siCtrl'],
                            #     y = figdata$rel_unpr[figdata$facetx == synchro & figdata$treatment == 'siAsf1']
                            # )

                            # kwres <- kruskal.test(rel_unpr ~ treatment, data = figdata[figdata$facetx == synchro,])
                            # print(kwres)
                        } else if (synchro == 'synchro') {
                            for (rel_time_v in c("0h", "2h", "5h")) {
                                message(paste0('Treatment effect on ', synchro, 'cells at ', rel_time_v, ' of Sphase:'))

                                glmres <- glm(
                                    data = figdata[figdata$facetx == synchro & figdata$rel_time == rel_time_v,],
                                    formula = rel_unpr ~ replicate + treatment + treatment:replicate
                                )
                                print(summary(glmres))

                                # kwres <- kruskal.test(rel_unpr ~ treatment, data = figdata[figdata$facetx == synchro & figdata$rel_time == rel_time_v,])
                                # print(kwres)
                            }
                        }

                    }
                }
            }

            ## the mean values for each boxplot (~condition)
            if (TRUE) {
                meanvals <- tapply(figdata$rel_unpr, figdata$condition, mean)
                meanvals[startsWith(names(meanvals), 'siAsf1')] / meanvals[startsWith(names(meanvals), 'siCtrl')]
                meanvals[grepl('_2h', names(meanvals))] / meanvals[grepl('_0h', names(meanvals))]
                meanvals[grepl('_5h', names(meanvals))] / meanvals[grepl('_2h', names(meanvals))]
            }
        }

        ## Ctrl vs siAsf1
        if (TRUE) {
            figdata <- mafdata[row_log, c('repl_histone', colnames(mafdata)[grepl('^si(Ctrl|Asf1)_', colnames(mafdata))])]
            figdata <- pivot_longer(figdata, colnames(figdata)[2:ncol(figdata)], names_to = 'condition', values_to = 'rel_unpr')
            figdata <- bind_cols(
                figdata,
                (tibble(data.frame(do.call(rbind, strsplit(figdata$condition, split='_')))) %>% mutate(treatment = X1, rel_time = X2))[c('treatment', 'rel_time')]
            )
            figdata$treatment <- factor(figdata$treatment, levels = c('siCtrl', 'siAsf1'), ordered = TRUE)
            figdata$facetx <- c('synchro', 'asyn', 'nasc')[(1 + (figdata$rel_time %in% c('asyn', 'nasc')) + (figdata$rel_time == 'nasc'))]

            ## plot Ctrl vs siAsf1 by synchronisation time
            for (syntime in names(ctrl_to_sirna_relation)) {
                if (syntime == "nasc") {
                    max_val <- max(as.matrix(mafdata[row_log, (ncol(mafdata) - 1):ncol(mafdata)])) * 1.1
                } else {
                    max_val <- max(as.matrix(mafdata[row_log, 7:(ncol(mafdata) - 2)])) * 1.1
                }
                xcol <- paste0("siCtrl_", syntime)
                ycol <- paste0("siAsf1_", syntime)

                figdata <- mafdata[row_log, c("repl_histone", xcol, ycol)]
                # figdata <- mafdata[row_log, c("repl_histone", paste0(xcol, "_", 1:2), paste0(ycol, "_", 1:2))]
                # figdata <- bind_rows(
                #     afdata[row_log, "repl_histone"],
                #     afdata[row_log, "repl_histone"]
                # )
                # figdata <- bind_cols(
                #     figdata,
                #     (afdata[row_log, paste0(xcol, "_", 1:2)] %>% pivot_longer(cols = paste0(xcol, "_", 1:2), names_to = "replicate", values_to = xcol))
                # )
                # figdata$replicate <- sub('siCtrl_', '', figdata$replicate)
                # figdata <- bind_cols(figdata, (afdata[row_log, paste0(ycol, "_", 1:2)] %>% pivot_longer(cols = paste0(ycol, "_", 1:2), names_to = "replicate", values_to = ycol))[ycol])

                # figdata[[xcol]] <- sqrt(figdata[[xcol]])
                # figdata[[ycol]] <- sqrt(figdata[[ycol]])

                print(cor.test(figdata[[ycol]], figdata[[xcol]]))
                regr_line <- lm(figdata[[ycol]]~figdata[[xcol]])
                print(regr_line)

                ## PDF
                if (TRUE) {
                    fig <- ggplot(
                        figdata,
                        aes_string(
                            x = xcol,
                            y = ycol
                        )
                    ) +
                        geom_point() +
                        geom_abline(intercept = 0, slope = 1) +
                        geom_abline(intercept = regr_line$coefficient[1], slope = regr_line$coefficient[2], col = "#AA0000") +
                        # scale_y_continuous(trans = 'sqrt') +
                        # xlim(0, max(figdata[[xcol]], figdata[[ycol]])) +
                        # ylim(0, max(figdata[[xcol]], figdata[[ycol]])) +
                        xlim(0, max_val) +
                        ylim(0, max_val) +
                        theme(text = element_text(size = 30))

                    figname <- paste0(figdir, '/', 'siCtrl_vs_siAsf1', '_', syntime, '_avg_fraction.pdf'); message(figname)
                    pdf(figname)
                    print(fig)
                    bouh <- dev.off()
                }

                ## HTML
                if (TRUE) {
                    hfigdata <- figdata
                    hfigdata[['siCtrl']] <- hfigdata[[xcol]]
                    hfigdata[['siAsf1']] <- hfigdata[[ycol]]

                    fig <- plot_ly(data = hfigdata,
                        type = 'scatter',
                        text = hfigdata$repl_histone
                    )
                    fig <- fig %>% add_trace(
                        x = ~siCtrl,
                        y = ~siAsf1,
                        # colors = rep(gg_color_hue(8), length(levels(resu$layout$organ)) %/% 8 + 1)[1:length(levels(resu$layout$organ))],
                        mode = 'markers'
                    ) #, size=0.1)

                    figname <- paste0(figdir, "/siCtrl_vs_siAsf1_", syntime, "_avg_fraction.html"); message(figname)
                    htmlwidgets::saveWidget(fig, figname, selfcontained = F)
                }
            }

            ## plot specifically for H1-2 gene
            if (TRUE) {
                tmpdata <- mafdata[mafdata$repl_histone == "H1-2",c('repl_histone', colnames(mafdata)[7:ncol(mafdata)])]
                figdata <- pivot_longer(tmpdata[c(1, 2:4, 8, 10)], colnames(tmpdata)[c(2:4, 8, 10)], names_to = "timing", values_to = "siCtrl")
                figdata[['siAsf1']] <- pivot_longer(tmpdata[c(1, 5:7, 9, 11)], colnames(tmpdata)[c(5:7, 9, 11)], names_to = "timing", values_to = "siAsf1")[["siAsf1"]]
                figdata$timing <- do.call(rbind, strsplit(figdata$timing, '_'))[,2]
                figdata$timing <- factor(figdata$timing, levels = figdata$timing)

                max_val <- max(unlist(tmpdata[2:ncol(tmpdata)]))

                fig <- ggplot(
                    figdata,
                    aes(
                        x = siCtrl,
                        y = siAsf1,
                        col = timing
                    )
                ) +
                    geom_point(size = 5) +
                    geom_abline(intercept = 0, slope = 1) +
                    geom_abline(mapping = aes(intercept = 0, slope = siAsf1 / siCtrl, col = timing)) +
                    xlim(0, max_val) +
                    ylim(0, max_val) +
                    theme(text = element_text(size = 30), axis.text.x = element_text(angle = 45, hjust = 1))

                figname <- paste0(figdir, '/', 'H1-2_avg_fraction_of_unclivedmRNA_2D.pdf'); message(figname)
                pdf(figname)
                print(fig)
                bouh <- dev.off()
            }

            ## Calculate statistical significance (by GLM on read counts)
            if (TRUE) {
                moddata <- fdata[row_log, c('repl_histone', colnames(fdata)[7:ncol(fdata)])] %>% pivot_longer(colnames(fdata)[7:ncol(fdata)], names_to = 'sample', values_to = 'count')
                moddata <- moddata[! moddata$sample %in% c(paste0('D356T', 17:20)),]
                moddata[['sqrt_count']] <- sqrt(moddata$count)
                
                moddata <- left_join(
                    moddata,
                    desdf,
                    by = 'sample'
                )
                # moddata <- moddata[! moddata$time_rel %in% c('0', 'asyn'),]
                # moddata <- moddata[moddata$time_rel %in% c("2", "5"),]

                if (FALSE) {
                    moddata <- afdata[row_log, c('repl_histone', colnames(afdata)[grepl('^siCtrl', colnames(afdata))])] %>% pivot_longer(colnames(afdata)[grepl('^siCtrl', colnames(afdata))], names_to = 'sample_siCtrl', values_to = 'count_siCtrl')
                    moddata <- tibble(cbind(moddata,
                        (afdata[row_log, c('repl_histone', colnames(afdata)[grepl('^siAsf1', colnames(afdata))])] %>% pivot_longer(colnames(afdata)[grepl('^siAsf1', colnames(afdata))], names_to = 'sample_siAsf1', values_to = 'count_siAsf1'))[,2:3]
                    ))
                    moddata[['slope']] <- (moddata$count_siAsf1 + 10^-5) / (moddata$count_siCtrl + 10^-5)
                    # hist(log10(moddata$slope))
                    moddata[['time_rel']] <- factor(sub('h$', '', do.call(rbind, strsplit(moddata$sample_siAsf1, '_'))[,2]), levels = c(0, 2, 5, "asyn"))
                }


                modelsimple <- sqrt_count ~ treatment
                glmres1 <- glm(
                    formula = modelsimple,
                    family = gaussian,
                    data = moddata
                )
                summary(glmres1)

                modelform <- sqrt_count ~ treatment + time_rel:treatment
                glmres2 <- glm(
                    formula = modelform,
                    family = gaussian,
                    data = moddata
                )
                summary(glmres2)
            }
        }

    }
}


####
