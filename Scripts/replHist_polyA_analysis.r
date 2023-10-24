
{
    require(tidyverse)
    require(rtracklayer)
    require(edgeR)
    require(data.table)
    require(reshape2)
}


#####################################################
#### Setting up environment
#####################################################
{
    ## general variables
    {
        oldtheme <- theme_set(theme_bw())
        theme_txt <- theme(text = element_text(size = 20))
        theme_txt_xangle <- theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust=1))

        datadir <- "./Data/ASF1_chaperone"
        resdir <- "./Results/ASF1_chaperone/polyA"; dir.create(resdir, showWarnings = F, recursive = F)
        figdir <- "./Images/ASF1_chaperone/polyA"; dir.create(figdir, showWarnings = F, recursive = F)

        uncliv_data_dir <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/RNA-Seq_fragment_on_HDE_SL_and_upstream'
        read_count_dir <- './Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/featureCounts/on_coding_sequence_of_replicative_histone_genes' #'/Volumes/Seb_GA/calcsub/slemaire/processed' #'./Data/ASF1_chaperone/RNA-Seq/read_counts'

        condition_colpal <- c('siCtrl' = '#3174a1', 'siAsf1' = '#c03d3e')

        gene_set_name_list <- list(
            'variant_histone_and_histone_chaperone_genes' = 'histVarChap',
            'replicative_histone_genes' = 'replHist',
            'decile_9' = 'dec9Expr'
        )
    }

    ## load the sets of genes
    {
        ## whole genome
        genes <- rtracklayer::import.gff(file.path('./Data/GRCh38.104.sorted.geneOnly.gtf'))

        ## histone & chaperone genes
        hc_tab <- tibble(read.table('./Data/GtexTcgaGeneExpr/histone_chaperone_complete_list.csv', h=T, sep = '\t'))
        hc_gff <- genes[genes$gene_id %in% hc_tab$EnsemblGeneId]

        ## only replicative histone genes
        rhc_tab <- hc_tab[hc_tab$subClass == 'replicative', ]
        rhc_gff <- hc_gff[hc_gff$gene_id %in% rhc_tab$EnsemblGeneId]

        ## get unique EnsEMBL annotation of coding sequenes for replicative histones genes
        cds_rhc_p <- './Data/ASF1_chaperone/annotations/replicative_histone_CDS.gtf'
        cdsgr <- rtracklayer::import.gff(cds_rhc_p)

        ## histone variants and histone chaperones
        hvc_tab <- hc_tab[hc_tab$subClass != 'replicative', ]
        hvc_gff <- hc_gff[hc_gff$gene_id %in% hvc_tab$EnsemblGeneId]

        ## gene of the second decile of most expressed genes (in synchronized cells)
        gdec_gff <- rtracklayer::import.gff("./Results/ASF1_chaperone/decile_expr_genes_increasing_syncCellsRNAseq/decile_9.gff")
    }

    ## Document duplicated replicative histones
    {
        dup.l <- list(
                c('H4C14', 'H4C15'), # c('HIST2H4A', 'HIST2H4B'),
                c('H3C14', 'H3C15'), # c('HIST2H3C', 'HIST2H3A'),
                c('H2AC18', 'H2AC19') # c('HIST2H2AA3', 'HIST2H2AA4')
        )

        dup.ll <- lapply(dup.l, function(xxx) {xxx[1]})
        names(dup.ll) <- lapply(dup.l, function(xxx) {xxx[2]})

        dupi.l <- lapply(dup.l, function(xxx) {
            left_join(
                data.frame('gene_name'=xxx),
                as.data.frame(mcols(genes))[c('gene_id', 'gene_name')],
                by='gene_name'
            )$gene_id
        })

        dupi.ll <- lapply(dupi.l, function(xxx) {xxx[1]})
        names(dupi.ll) <- lapply(dupi.l, function(xxx) {xxx[2]})
    }

    ## load counts of poly-A read pairs
    {
        polya_p <- "./Results/ASF1_chaperone/poly_adenylation_22T/polyA_22T_readPair_counts.txt"
        pa_tab <- tibble(read.table(polya_p, h=F))
        colnames(pa_tab) <- c('sample', 'gene_set', 'count')

        pa_tab$sample <- factor(pa_tab$sample, levels=unique(pa_tab$sample))

        pa_tab$gene_set_lab <- as.vector(unlist(gene_set_name_list[pa_tab$gene_set]))
    }

    ## RNA-seq design
    {
        spll <- list(
            "D356T1" = "siCtrl_0h_R1",
            "D356T2" = "siCtrl_0h_R2",
            "D356T3" = "siCtrl_2h_R1",
            "D356T4" = "siCtrl_2h_R2",
            "D356T5" = "siCtrl_5h_R1",
            "D356T6" = "siCtrl_5h_R2",
            "D356T7" = "siAsf1_0h_R1",
            "D356T8" = "siAsf1_0h_R2",
            "D356T9" = "siAsf1_2h_R1",
            "D356T10" = "siAsf1_2h_R1",
            "D356T11" = "siAsf1_5h_R1",
            "D356T12" = "siAsf1_5h_R2",
            "D356T13" = "siCtrl_asyn_total_RNA_R1",
            "D356T14" = "siCtrl_asyn_total_RNA_R2",
            "D356T15" = "siAsf1_asyn_total_RNA_R1",
            "D356T16" = "siAsf1_asyn_total_RNA_R2",
            "D356T17" = "siCtrl_asyn_nascent_RNA_R1",
            "D356T18" = "siAsf1_asyn_nascent_RNA_R1",
            "D356T19" = "siCtrl_asyn_nascent_RNA_R2",
            "D356T20" = "siAsf1_asyn_nascent_RNA_R2"
        )

        ## experimental design
        desdf <- tibble(data.frame(
            "sample" = names(spll),
            "condition" = factor(c(rep("siCtrl", 6), rep("siAsf1", 6), rep("siCtrl", 2), rep("siAsf1", 2), "siCtrl", "siAsf1", "siCtrl", "siAsf1"), levels = c("siCtrl", "siAsf1"), order = F),
            "Sphase" = factor(c(rep(c(0, 0, 2, 2, 5, 5), 2), rep("asyn", 8)), levels = c(0, 2, 5), order = F),
            "RNA" = factor(c(rep("total", 16), rep("nascent", 4)), levels = c("total", "nascent"), order = F),
            "replicate" = c(rep(1:2, 8), rep(1:2, each=2))
        ))

        desdf$spl_name <- as.vector(unlist(spll[desdf$sample]))
        desdf$setting <- sub('(_RNA|)_R[0-9]$', '', desdf$spl_name)
        desdf$setting <- factor(desdf$setting, levels=c("siCtrl_0h", "siCtrl_2h", "siCtrl_5h", "siAsf1_0h", "siAsf1_2h", "siAsf1_5h", "siCtrl_asyn_total", "siAsf1_asyn_total", "siCtrl_asyn_nascent", "siAsf1_asyn_nascent"))
    }

    ## load spike-in x sequencing depth normalization factors
    if (TRUE) {
        spi.seqdep.nf.t <- read.table(file.path('./Results/ASF1_chaperone', 'spike_and_seqDepth_norm_factor.tsv'), h=T, sep='\t')
        spi.seqdep.nf <- as.vector(spi.seqdep.nf.t$norm_factor)
        names(spi.seqdep.nf) <- spi.seqdep.nf.t$sample
    }

    ## Load the clusters by Alberto
    {
        albcl.p <- './Data/ASF1_chaperone/Alberto_Gatto_results/results/clusters_6/gene_clusters_sync.tsv'
        albcl.t <- read.table(albcl.p, h=T, sep='\t')
        colnames(albcl.t)[1] <- 'gene_id'
        # head(albcl.t)

        ## combine the duplicated replicative histone genes in the clustering
        for (dgi in 1:length(dupi.l)) {
            # dgi <- 2
            dgp <- dupi.l[[dgi]]
            
            salbcl.t <- albcl.t[albcl.t$gene_id %in% dgp,]
            # print(stab)

            if (2 %in% salbcl.t$cluster) {
                albcl.t$cluster[albcl.t$gene_id %in% dgp] <- 2
                albcl.t <- albcl.t[!(albcl.t$gene_id %in% dgp[2]),]
                # str(albcl.t)
            }
        }
    }
}


#####################################################
#### Analyze contingencies of "polyA" reads
#####################################################
{
    ## load counts of gene expression
    {
        ## compute path to read count files
        rcdir <- "./Data/ASF1_chaperone/RNA-Seq/read_counts"
        splpl <- paste0(rcdir, '/', names(spll), '_geneByExon.rv')

        ## load the read count tables
        cntl <- list()
        for (iii in 1:length(splpl)) {
            spl <- names(spll)[iii]
            cntl[[spl]] <- fread(splpl[[iii]], h = T)[,c(1, 7)]
            colnames(cntl[[spl]]) <- c("Geneid", "count")
        }

        ## merge read counts in one file
        mcnt <- cntl[[1]][,"Geneid"]
        for (spl in names(cntl)) {
            mcnt <- cbind(mcnt, cntl[[spl]][,"count"])
        }
        colnames(mcnt) <- c("Geneid", names(cntl))

        mcnt <- tibble(mcnt)
        cnt_col <- which(grepl('^D356T', colnames(mcnt)))

        ## combine the duplicated replicative histone genes
        for (dgi in 1:length(dupi.l)) {
            # dgi <- 1
            dgp <- dupi.l[[dgi]]
            
            smcnt <- mcnt[mcnt$Geneid %in% dgp,]
            # print(stab)

            sel_col <- !grepl('^(Geneid)', colnames(smcnt))
            smcnt[1, sel_col] <- smcnt[1, sel_col] + smcnt[2, sel_col]

            mcnt[mcnt$Geneid %in% dgp[1], sel_col] <- smcnt[1, sel_col]
            mcnt <- mcnt[!(mcnt$Geneid %in% dgp[2]),]
        }

        ## density of genes per total read count
        if (FALSE) {
            ggplot(
                data.frame('total_read_count' = rowSums(mcnt[cnt_col])),
                aes(
                    x = total_read_count + 1
                )
            ) +
                scale_x_continuous(trans='log10') +
                geom_density()
        }


        ## filter out low saturation genes
        dge <- DGEList(as.matrix(mcnt[cnt_col]))
        keep <- filterByExpr(dge)
        filt.mcnt <- mcnt[keep,]

        ## normalize by gene length
        filt.mcnt$gene_length <- left_join(
            filt.mcnt['Geneid'],
            tibble(as.data.frame(genes)) %>% mutate(Geneid=gene_id),
            by = 'Geneid'
        )[['width']]

        gl.mcnt <- filt.mcnt
        for (spl_name in colnames(filt.mcnt)[cnt_col]) {
            gl.mcnt[[spl_name]] <- gl.mcnt[[spl_name]] / gl.mcnt$gene_length
        }

        ## normalize by TMM algorithm
        dge <- DGEList(as.matrix(gl.mcnt[cnt_col]))
        dge <- calcNormFactors(dge)
        rawCPM <- cpm(dge, log=FALSE)
        # logCPM <- cpm(dge, log=TRUE)

        tmm.gl.mcnt <- gl.mcnt
        for (spl_name in colnames(mcnt)[cnt_col]) {
            tmm.gl.mcnt[[spl_name]] <- rawCPM[, spl_name]
        }

        ## mean vs sd of expression across samples
        if (FALSE)
        {
            md.tab <- tibble(
                data.frame(
                    'mean_expr' = rowMeans(tmm.gl.mcnt[cnt_col]),
                    'std_expr' = apply(tmm.gl.mcnt[cnt_col], 1, sd)
                )
            )

            ggplot(
                md.tab,
                aes(
                    x = mean_expr,
                    y = std_expr
                )
            ) +
                scale_x_continuous(trans="log10") + 
                scale_y_continuous(trans="log10") + 
                geom_point()

            ggplot(
                md.tab,
                aes(
                    x = mean_expr
                )
            ) +
                scale_x_continuous(trans="log10") + 
                geom_density()
        }

    }


    ## Rate of poly-A RNAs in total expression (by transcript point of view)
    {
        ## link gene set with label
        gene_set_list <- list(
            'histVarChap' = hvc_gff$gene_id,
            'replHist' = rhc_gff$gene_id,
            'dec9Expr' = gdec_gff$gene_id
        )

        ## recover in each sample the total expression per gene set
        totExprGeneSet_list <- list()
        for (sample in desdf$sample) {
            totExprGeneSet_list[[sample]] <- list()
            for (gene_set in pa_tab$gene_set_lab) {
                totExprGeneSet_list[[sample]][[gene_set]] <- sum(tmm.gl.mcnt[[sample]][tmm.gl.mcnt$Geneid %in% gene_set_list[[gene_set]]])
            }
        }

        totExprGeneSet <- tibble(reshape2::melt(totExprGeneSet_list))
        colnames(totExprGeneSet) <- c('expr', 'gene_set_lab', 'sample')

        ## normalize the poly-A counts
        {
            totExprGeneSet_mat <- matrix(totExprGeneSet$expr, nrow=3, dimnames=list(unique(totExprGeneSet$gene_set_lab), unique(totExprGeneSet$sample)))
            polyA_mat <- matrix(pa_tab$count, nrow=3, dimnames=list(unique(pa_tab$gene_set_lab), unique(pa_tab$sample)))

            norm_polyA_mat <- polyA_mat / totExprGeneSet_mat

            pa_tab$norm_count <- as.vector(norm_polyA_mat) * 1000
        }

        ## plot the normalized poly-A conuts
        {
            gdata <- pa_tab
            gdata <- left_join(
                gdata,
                desdf[c('sample', 'setting', 'condition')],
                by = 'sample'
            )

            fig <- ggplot(
                gdata,
                aes(
                    x = setting,
                    y = norm_count,
                    fill = condition,
                    col = "black",
                    group = sample
                )
            ) +
                facet_wrap(~gene_set_lab, scale = 'free_y') +
                geom_col(position = "dodge") +
                # geom_point() +
                theme_txt_xangle +
                scale_color_manual(values=c("black"="#000000")) +
                scale_fill_manual(values=condition_colpal)
            
            figname <- file.path(figdir, 'polyA_rate_in_total_expression_barplot.svg'); message(figname)
            svg(figname, width = 3 * 7)
            print(fig)
            bouh <- dev.off()
        }
    }


    ## Rate of poly-A RNAs in total expression (by reads in gene annotation point of view)
    {
        ## Compile the counts of reads per replicative histone
        rtab <- list()
        seq_depth_list <- list()
        for (filename in list.files(read_count_dir)) {
            if (grepl('replicative_histone_CDS_D356T[0-9]*.tsv$', filename)) {
                message(filename)
                intab <- read.table(paste0(read_count_dir, '/', filename), h = T)
                seq_depth_list[[filename]] <- sum(intab[[7]])
                intab <- intab[intab$Geneid %in% cdsgr$gene_id,]
                rtab[[filename]] <- intab
            }
        }

        ## Compile the total mRNA read counts for the replication histone genes
        {
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

            ## combine the duplicated replicative histone genes
            for (dgi in 1:length(dupi.l)) {
                # dgi <- 1
                dgp <- dupi.l[[dgi]]
                
                srcdata <- rcdata[rcdata$Geneid %in% dgp,]
                # print(stab)

                sel_col <- !grepl('^(Geneid|Chr|Start|End|Strand|Length|repl_histone)', colnames(srcdata))
                srcdata[1, sel_col] <- srcdata[1, sel_col] + srcdata[2, sel_col]

                rcdata[rcdata$Geneid %in% dgp[1], sel_col] <- srcdata[1, sel_col]
                rcdata <- rcdata[!(rcdata$Geneid %in% dgp[2]),]
            }

            ## normalize the counts by annotation length
            rcdata_lenNorm <- rcdata
            for (coln in colnames(rcdata_lenNorm)[grepl('^D356T', colnames(rcdata_lenNorm))]) {
                rcdata_lenNorm[[coln]] <- rcdata_lenNorm[[coln]] / rcdata_lenNorm$Length
            }
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
            npa_tab <- pa_tab
            nrcdata_lenNorm <- rcdata_lenNorm

            for (coln in cat.dds.ld$sample) {
                npa_tab[['norm_count']][npa_tab$sample == coln] <- npa_tab[['count']][npa_tab$sample == coln] / cat.dds.ld$sizeFactor[cat.dds.ld$sample == coln]
                nrcdata_lenNorm[[coln]] <- nrcdata_lenNorm[[coln]] / cat.dds.ld$sizeFactor[cat.dds.ld$sample == coln]
            }

            for (coln in nas.dds.ld$sample) {
                npa_tab[['norm_count']][npa_tab$sample == coln] <- npa_tab[['count']][npa_tab$sample == coln] / nas.dds.ld$sizeFactor[nas.dds.ld$sample == coln]
                nrcdata_lenNorm[[coln]] <- nrcdata_lenNorm[[coln]] / nas.dds.ld$sizeFactor[nas.dds.ld$sample == coln]
            }

            nrcdata_lenNorm <- nrcdata_lenNorm[7:ncol(nrcdata_lenNorm)]
        }

        ## normalize by the number of unclived RNAs
        {
            repl_npa_tab <- npa_tab[npa_tab$gene_set_lab == 'replHist',]

            tot_cds <- colSums(nrcdata_lenNorm[repl_npa_tab$sample])
            repl_npa_tab$tot_cds <- as.vector(tot_cds)

            repl_npa_tab$norm_count_in_cds <- repl_npa_tab$count / repl_npa_tab$tot_cds
        }

        ## plot the rate of poly-A in unclived RNAs
        {
            gdata <- repl_npa_tab
            gdata <- left_join(
                gdata,
                desdf[c('sample', 'setting', 'condition')],
                by = 'sample'
            )

            fig <- ggplot(
                gdata,
                aes(
                    x = setting,
                    y = norm_count_in_cds,
                    fill = condition,
                    col = "black",
                    group = sample
                )
            ) +
                # facet_wrap(~gene_set_lab, scale = 'free_y') +
                geom_col(position = "dodge") +
                # geom_point() +
                theme_txt_xangle +
                scale_color_manual(values=c("black"="#000000")) +
                scale_fill_manual(values=condition_colpal)
            
            figname <- file.path(figdir, 'polyA_rate_in_CDS_barplot.svg'); message(figname)
            svg(figname, width = 1 * 7)
            print(fig)
            bouh <- dev.off()
        }
    }


    ## Rate of poly-A RNAs in unclived RNAs (can only be applied to replicative histoness)
    {
        ## Load the counts of unclived transcripts (~fragments)
        {
            ctab <- list()
            for (filename in list.files(uncliv_data_dir)) {
                if (grepl('D356T[0-9]*\\.sorted\\.markdup\\.mapq10\\.Nfiltered\\.hSL\\.20upstrOfSL\\.HDE\\.frag500bMax\\.countPerHDE\\.bed\\.gz', filename)) {
                    message(filename)
                    intab <- read.table(paste0(uncliv_data_dir, '/', filename))
                    colnames(intab) <- c('ref', 'pos_0.based', 'end', 'trash', 'repl_histone', 'strand', 'count')
                    ctab[[filename]] <- intab
                }
            }

            unc_tab <- tibble(ctab[[1]])
            unc_tab <- unc_tab[1:ncol(unc_tab) - 1]

            for (name in names(ctab)) {
                unc_tab[[strsplit(name, '\\.')[[1]][1]]] <- ctab[[name]][[7]]
            }

            unc_cnt_col <- which(grepl('^D356T', colnames(unc_tab)))
        }

        ## combine the duplicated replicative histone genes
        for (dgi in 1:length(dupi.l)) {
            # dgi <- 1
            dgp <- dup.l[[dgi]]
            
            sunc_tab <- unc_tab[unc_tab$repl_histone %in% dgp,]
            # print(stab)

            sel_col <- !grepl('^(ref|pos_0.based|end|trash|repl_histone|strand)', colnames(sunc_tab))
            sunc_tab[1, sel_col] <- sunc_tab[1, sel_col] + sunc_tab[2, sel_col]

            unc_tab[unc_tab$repl_histone %in% dgp[1], sel_col] <- sunc_tab[1, sel_col]
            unc_tab <- unc_tab[!(unc_tab$repl_histone %in% dgp[2]),]
        }

        ## normalize by the number of unclived RNAs
        {
            repl_pa_tab <- pa_tab[pa_tab$gene_set_lab == 'replHist',]

            tot_uncl <- colSums(unc_tab[repl_pa_tab$sample])
            repl_pa_tab$unclived <- as.vector(tot_uncl)

            repl_pa_tab$norm_count_in_unc <- repl_pa_tab$count / repl_pa_tab$unclived
        }

        ## plot the rate of poly-A in unclived RNAs
        {
            gdata <- repl_pa_tab
            gdata <- left_join(
                gdata,
                desdf[c('sample', 'setting', 'condition')],
                by = 'sample'
            )

            fig <- ggplot(
                gdata,
                aes(
                    x = setting,
                    y = norm_count_in_unc,
                    fill = condition,
                    col = "black",
                    group = sample
                )
            ) +
                # facet_wrap(~gene_set_lab, scale = 'free_y') +
                geom_col(position = "dodge") +
                # geom_point() +
                theme_txt_xangle +
                scale_color_manual(values=c("black"="#000000")) +
                scale_fill_manual(values=condition_colpal)
            
            figname <- file.path(figdir, 'polyA_rate_in_unclived_barplot.svg'); message(figname)
            svg(figname, width = 1 * 7)
            print(fig)
            bouh <- dev.off()
        }
    }
}


#####################################################
#### Analyze contingencies of PAS (AAUAAA) in "polyA" reads
#####################################################

resdir_polya <- './Results/ASF1_chaperone/poly_adenylation_22T'

## Evaluate the absolute occurence of PAS in "polyA" reads
samt.ll <- list()
samgi.ll <- list()
has.pas.nbr.ll <- list()
for (gene_set in c('replicative_histone_genes', 'variant_histone_and_histone_chaperone_genes', 'decile_9'))
{
    # gene_set <- 'replicative_histone_genes'
    
    ## load the identified polyA-reads
    {
        samt.l <- list()
        samgi.l <- list()
        for (spln in desdf$sample) {
            samf <- file.path(resdir_polya, paste0(gene_set, '_', spln, '_polyAreads.sam'))
            message(samf)
            tryCatch({
                samt <- read.table(samf, h=F, sep='\n', comment.char = "@")
                samt <- as.list(samt[[1]])
            }, error = function(e) {
                samt <- list()
            })

            samt <- lapply(samt, function(xxx) {
                strsplit(xxx, split='\t')[[1]]
            })

            ## filter out replicative histone genes not in group 2
            if (gene_set == 'replicative_histone_genes') {
                gene.id.vec <- unlist(lapply(samt, function(xxx) {
                    strsplit(tail(xxx, n=1), split='GI:')[[1]][2]
                }))

                samt <- samt[!c(gene.id.vec %in% c('HIST2H4A', 'HIST2H4B', 'HIST1H2AA', 'HIST1H2BA', 'HIST1H1A', 'HIST1H2BB', 'HIST1H3C', 'HIST1H1T', 'HIST1H2BE', 'HIST1H3E', 'HIST1H2BM', 'HIST1H3I'))]
                # samt <- samt[gene.id.vec %in% albcl.t$gene_id[albcl.t$cluster == 2]]
            }

            samgi <- lapply(samt, function(xxx) {
                yyy <- tail(xxx, n=1)
                yyy <- tail(strsplit(yyy, split=' ')[[1]], n=1)
                yyy <- strsplit(yyy, split=':')[[1]][2]
            })
            samgi <- unlist(samgi)
            samgi[samgi %in% names(dupi.ll)] <- unlist(dupi.ll[samgi[samgi %in% names(dupi.ll)]])

            samt.l[[spln]] <- samt
            samgi.l[[spln]] <- samgi
        }

        samt.ll[[gene_set_name_list[[gene_set]]]] <- samt.l
        samgi.ll[[gene_set_name_list[[gene_set]]]] <- samgi.l
    }


    ## look for PAS 'AAUAAA' in the polyA read sequences
    {
        has.pas.nbr.l <- list()
        for (spln in names(samt.l)) {
            samt <- samt.l[[spln]]
            
            has.pas <- lapply(samt, function(xxx) {
                # grepl('TTTATT', xxx[10])
                TRUE
            })

            has.pas.nbr <- sum(unlist(has.pas))
            has.pas.nbr.l[[spln]] <- has.pas.nbr
        }
        has.pas.nbr.l <- unlist(has.pas.nbr.l)
        has.pas.nbr.l <- data.frame(
            'sample' = names(has.pas.nbr.l),
            'count' = as.vector(has.pas.nbr.l)
        )

        has.pas.nbr.ll[[gene_set_name_list[[gene_set]]]] <- has.pas.nbr.l
    }


    ## build plot of PAS contingencies
    {
        # per dataset and condition

        gdata <- has.pas.nbr.l
        gdata <- left_join(
            gdata,
            desdf,
            by='sample'
        )
        gdata$group <- paste0(gdata$RNA, '_', gdata$Sphase, '_', gdata$condition)
        gdata$group <- factor(gdata$group, levels=unique(gdata$group))
        gdata$experiment <- 'Sphase'
        gdata$experiment[is.na(gdata$Sphase)] <- as.character(gdata$RNA[is.na(gdata$Sphase)])
        gdata$experiment <- factor(gdata$experiment, levels=unique(gdata$experiment))

        seg.gdata <- pivot_wider(gdata %>% mutate(repl_names = paste0('repl', replicate)), id_cols=c('group', 'experiment', 'condition'), names_from = 'repl_names', values_from = 'count')
        seg.gdata$avg <- apply(seg.gdata[c('repl1', 'repl2')], 1, mean)

        fig <- ggplot(
            gdata,
            aes(
                x = group,
                y = count,
                group = group,
                fill = condition
            )
        ) +
            facet_wrap(~experiment, scales = "free_x") +
            geom_col(
                data = seg.gdata,
                mapping = aes(
                    y = avg,
                    fill = condition
                ),
            ) +
            geom_segment(
                data = seg.gdata,
                mapping = aes(
                    y = repl1,
                    xend = group,
                    yend = repl2
                )
            ) +
            geom_point() + 
            scale_fill_manual(values=condition_colpal) +
            theme_txt_xangle
    }


    ## print the plot
    {
        figname <- file.path(figdir, paste0('AAUAAA_contingency_in_polyAreads_', gene_set, '.pdf')); message(figname)
        pdf(figname)
        print(fig)
        bouh <- dev.off()
    }
}


## Evaluate the fraction of PAS in "polyA" reads normalized by spike-in
for (gene_set in c('replicative_histone_genes', 'variant_histone_and_histone_chaperone_genes', 'decile_9'))
{
    # gene_set <- 'replicative_histone_genes'

    ## compute fractions of PAS
    {
        ## combine total expression counts and PAS counts
        {
            has.pas.nbr.l <- has.pas.nbr.ll[[gene_set_name_list[[gene_set]]]]

            has.pas.fra.l <- has.pas.nbr.l
            has.pas.fra.l$normFact <- spi.seqdep.nf[has.pas.fra.l$sample]
            has.pas.fra.l$fraction <- has.pas.fra.l$count * has.pas.fra.l$normFact
        }
    }

    ## prepare data
    {
        gdata <- has.pas.fra.l
        gdata <- left_join(
            gdata,
            desdf,
            by='sample'
        )
        gdata$group <- paste0(gdata$RNA, '_', gdata$Sphase, '_', gdata$condition)
        gdata$group <- factor(gdata$group, levels=unique(gdata$group))
        gdata$experiment <- 'Sphase'
        gdata$experiment[is.na(gdata$Sphase)] <- as.character(gdata$RNA[is.na(gdata$Sphase)])
        gdata$experiment <- factor(gdata$experiment, levels=unique(gdata$experiment))

        seg.gdata <- pivot_wider(gdata %>% mutate(repl_names = paste0('repl', replicate)), id_cols=c('group', 'experiment', 'condition'), names_from = 'repl_names', values_from = 'fraction')
        seg.gdata$avg <- apply(seg.gdata[c('repl1', 'repl2')], 1, mean)

        log_nasc_spl <- gdata$sample %in% desdf$sample[desdf$RNA == 'nascent']
    }

    ## Plot of PAS contingencies
    {
        ## build figure for total RNA-seq datasets
        {
            sgdata <- gdata[!log_nasc_spl,]
            seg.sgdata <- seg.gdata[seg.gdata$experiment != 'nascent',]

            fig <- ggplot(
                sgdata,
                aes(
                    x = group,
                    y = fraction,
                    group = group,
                    fill = condition
                )
            ) +
                facet_wrap(~experiment, scales = "free_x") +
                geom_col(
                    col='#000000',
                    position = position_dodge2(width = 0.35, preserve = "single")
                ) +
                scale_fill_manual(values=condition_colpal) +
                theme_txt_xangle
        }

        ## print the figure for total RNA-seq datasets
        {
            figname <- file.path(figdir, paste0('AAUAAA_fractionInTotalExpr_in_polyAreads_', gene_set, '_totalRNAdatasets.pdf')); message(figname)
            pdf(figname)
            print(fig)
            bouh <- dev.off()
        }

        ## build figure for nascent RNA
        {
            sgdata <- gdata[log_nasc_spl,]
            seg.sgdata <- seg.gdata[seg.gdata$experiment == 'nascent',]

            fig <- ggplot(
                sgdata,
                aes(
                    x = group,
                    y = fraction,
                    group = group,
                    fill = condition
                )
            ) +
                facet_wrap(~experiment, scales = "free_x") +
                geom_col(
                    col='#000000',
                    position = position_dodge2(width = 0.35, preserve = "single")
                ) +
                scale_fill_manual(values=condition_colpal) +
                theme_txt_xangle
        }

        ## print the figure for nascent RNA
        {
            figname <- file.path(figdir, paste0('AAUAAA_fractionInTotalExpr_in_polyAreads_', gene_set, '_nascentRNA.pdf')); message(figname)
            pdf(figname, width = 0.5 * 7)
            print(fig)
            bouh <- dev.off()
        }
    }

}



####

named_vec_to_df <- function(nvec, id_name=NULL, value_name=NULL) {
    ## convert a vector of named item into a dataframe
    out <- tibble(
        data.frame(
            "name" = names(nvec),
            "value" = as.vector(nvec)
        )
    )

    if (!is.null(id_name)) {
        colnames(out)[1] <- id_name
    }

    if (!is.null(value_name)) {
        colnames(out)[2] <- value_name
    }

    return(out)
}


xxx <- samt.ll[["replHist"]] #['D356T14']

yyy <- lapply(samt.ll[["replHist"]], function(xxx) {

    gil <- lapply(xxx, function(zzz){
        if (grepl('TTTATT', zzz[10])) {
            gene_id <- strsplit(tail(zzz, n=1), split='GI:')[[1]][2]
            return(gene_id)
        }
    })
    gil <- unlist(gil)

    out <- named_vec_to_df(table(gil), id_name='gene_id')
    out <- left_join(
        out,
        albcl.t[c('gene_id', 'name')],
        by='gene_id'
    )    
    return(out)
})


####
