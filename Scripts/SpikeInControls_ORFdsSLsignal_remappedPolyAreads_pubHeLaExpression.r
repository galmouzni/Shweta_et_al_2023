## Differential gene expression analysis of Schweta's project (+/- ASF1)

{
    ## Should we update the recorded results ('save_data_f' functions)
    SAVE_DATA_LOG <- FALSE
    
    ## Shoul we print the plots
    PRINT_PLOTS <- FALSE
}


## PACKAGES
if (TRUE) {
    require(reshape2)
    require(data.table)
    require(tidyverse)
    require(ggdendro)
    require(grid)
    require(edgeR)
    require(DESeq2)
    library("vsn")
    library("pheatmap")
    library("RColorBrewer")
    library("IHW")
    require(igraph)
    require(WebGestaltR)
    library("biomaRt")
}

oldtheme <- theme_set(theme_bw())
theme_txt <- theme(text = element_text(size = 20))
theme_txt_xangle <- theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))
theme_txt_xvert <- theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1))


## FUNCTIONS
if (TRUE) {
    source('./Scripts/general/_functions_in_R.r')

    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }


    dds.analysis <- function(mcnt, desdf, design, Geneid.vec, thres.ind.filter = 10) {
        ## DESeq2 setting
        dds <- DESeqDataSetFromMatrix(countData = mcnt,
                                    colData = desdf,
                                    design = design)


        mcols(dds) <- DataFrame(mcols(dds), data.frame(Geneid = Geneid.vec))

        dds$condition <- factor(dds$condition, levels = c("siCtrl", "siAsf1"))

        ## independent filtering
        keep <- rowSums(counts(dds)) >= thres.ind.filter
        dds <- dds[keep,]


        ## DESeq processing
        dds <- DESeq(dds)

        ## transforming data (VST)
        vsd <- vst(dds, blind=FALSE)

        return(list("dds"=dds, "vsd"=vsd))
    }


    deg.summary.fc <- function(resl, expr.mtx, fdr.thresh = 0.01, scaled=TRUE) {
        ## adjust notation (Note: VST is not applied on the expression matrix in this function)
        vsd <- expr.mtx

        ##  # DEGs gene idx, up/down separated
        degil <- lapply(resl, function(xxx) {
            list(
                "up" = (!is.na(xxx$padj)) & xxx$padj <= 0.01 & xxx$log2FoldChange > 0,
                "down" = (!is.na(xxx$padj)) & xxx$padj <= 0.01 & xxx$log2FoldChange < 0
                )
        })

        ##  # merge DEG idx in one vector
        mdegil <- lapply(degil, function(xxx) {
            xxx[["up"]] | xxx[["down"]]
        })

        for (iii in 1:length(mdegil)) {
            mdegil[[1]] <- mdegil[[1]] | mdegil[[iii]]
        }
        mdegil <- mdegil[[1]]

        ##  # identify the tendencies (hierarchical clustering-ordered heatmap display)
        degvsd <- vsd[mdegil]
        degvsd.center <- t(apply(assay(degvsd), 1, ifelse(scaled,
            function(xxx) {
                (xxx - mean(xxx)) / sd(xxx)
            },
            function(xxx) {
                (xxx - mean(xxx))
            }
        )))

        deg.hclust <- hclust(d = dist(x = (as.matrix(degvsd.center))))
        deg.dendro <- as.dendrogram(deg.hclust)

        return(list(
            "res.deg.only"=degil,
            "deg.log"=mdegil,
            "deg.expr.mtx"=degvsd,
            "center.expr.mtx"=degvsd.center,
            "hclust"=deg.hclust,
            "dendrogram"=deg.dendro
        ))
    }


    ora_go_enrichment <- function(deg.hclust, degids, allids, nclust = 15, gedir = './geneEnrichAnalysis', enrichDatabaseVec=NA) {
        ## give default databases if not specified
        if (is.na(enrichDatabaseVec)) {
            enrichDatabaseVec <- c('geneontology_Biological_Process_noRedundant', 'geneontology_Cellular_Component_noRedundant', 'pathway_KEGG', 'pathway_Wikipathway')
        }

        ##
        deg.dendro <- as.dendrogram(deg.hclust)

        ## Cluster the DEG by their expression pattern
        if (TRUE) { 
            deg.grp <- cutree(tree = deg.hclust, k = nclust)
            grp.cnt <- table(deg.grp)
            print(grp.cnt)

            ## enrichment analysis by DEG cluster
            if (TRUE) {
                dir.create(gedir, showWarnings = F)
                erchl <- list()
                for (cii in unique(c(deg.grp[order.dendrogram(deg.dendro)]))) {
                    # cii <- 1
                    print(cii)
                    if (grp.cnt[as.character(cii)] < 10) {
                        message('Too few genes. Pass.')
                        next
                    }

                    for (enrdat in enrichDatabaseVec) {
                        ## ORA (for enrichment in a set of genes compared to a full set)
                        erchl[[as.character(cii)]][[enrdat]] <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                            enrichDatabase=enrdat,
                            interestGene=degids[deg.grp == cii],
                            interestGeneType="ensembl_gene_id",
                            referenceGene=mcols(vsd)$Geneid,
                            referenceGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=100, minNum=5,
                            outputDirectory=gedir
                        )
                    }
                }
            }
        }

        return(list(
            "deg.cluster.id"=deg.grp,
            "go.enrichment"=erchl
        ))
    }


    gsea_go_enrichment <- function(dds, fdr.thresh=0.05, enrichDatabaseVec=NULL, outputDirectory='.', ...) {
        if (is.null(enrichDatabaseVec)) {
            enrichDatabaseVec <- c('geneontology_Biological_Process_noRedundant', 'geneontology_Cellular_Component_noRedundant', 'pathway_KEGG', 'pathway_Wikipathway')
        }

        # get results from DESeq2 analysis
        res <- results(dds, ...)

        # get DEGs
        deg.log <- (!is.na(res$padj)) & (res$padj < fdr.thresh)

        # filter results data.frame for DEGs
        deg.res <- res[deg.log, ]

        # filter gene ids for DEGs
        deg.ids <- mcols(dds)$Geneid[deg.log]
        
        # GO enrichment
        erchl <- list()
        for (enrdat in enrichDatabaseVec) {
            # enrdat <- enrichDatabaseVec[1]

            erchl[[enrdat]] <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                enrichDatabase=enrdat,
                interestGene=data.frame("gene"=deg.ids, "score"=deg.res$log2FoldChange),
                interestGeneType="ensembl_gene_id",
                sigMethod="top", topThr=100, minNum=5,
                outputDirectory='.'
            )
        }

        return(erchl)
    }


    go_term_enrich_genes <- function(erchl, erch.fdr = 0.05, go.db = "geneontology_Biological_Process_noRedundant", mart = mart) {
        ## Retrieve genes contributing to the significant GO terms (Note: each gene is associated with the most significant (p-value) GO term it is associated with)
        resdf <- list()
        for (clustid in names(erchl)) {
            # clustid <- "14"
            erch <- erchl[[clustid]][[go.db]]

            go.sig.log <- erch$FDR < erch.fdr
            if (sum(go.sig.log)) {
                go.sig <- erch$geneSet[go.sig.log]

                resl <- lapply(as.list(erch$userId[go.sig.log]), function(eidv) {
                    ovlg <- strsplit(eidv, ';')[[1]]
                    results <- tibble(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                                    filters = "ensembl_gene_id", values = ovlg,
                                    mart = mart))
                })
                names(resl) <- erch$description[go.sig.log]

                # combine all GO term-associated genes (ordered by significativity of the GO terms)
                resdfe <- tibble(reshape2::melt(resl, id.vars=colnames(resl[[1]])))
                colnames(resdfe) <- c(head(colnames(resdfe), -1), "GO_term")

                # remove duplicated genes (keeping then associated with the most significant GO term)
                resdfe <- resdfe[!duplicated(resdfe$ensembl_gene_id), ]

                resdf[[clustid]] <- resdfe
            }
        }

        resdf <- tibble(reshape2::melt(resdf, id.vars = colnames(resdf[[1]])))
        colnames(resdf) <- c(head(colnames(resdf), -1), "cluster_id")

        return(resdf)
    }


    volcano_plots <- function(resl, volcdir) {
        for (compn in names(resl)) {
            # compn <- names(cat.resl)[1]
            resob <- resl[[compn]]

            maxlfc <- max(abs(resob$log2FoldChange[!is.na(resob$padj)]), na.rm = T)
            
            fig <- ggplot(
                tibble(data.frame(resob)) %>% filter(
                    !is.na(padj)
                ) %>% mutate(
                    sig = ifelse(padj < 0.05,
                        ifelse(
                            padj < 0.01,
                            ifelse(
                                padj < 0.001,
                                "< 0.001",
                                "< 0.01"
                            ),
                            "< 0.05"
                        ),
                        ">= 0.05"
                    )
                ),
                aes(
                    x = log2FoldChange,
                    y = padj,
                    col = sig
                )
            ) +
                scale_y_continuous(trans = "log10") +
                geom_point(size = 0.01) +
                xlim(- maxlfc, maxlfc) +
                theme(text = element_text(size = 30)) +
                guides(colour=guide_legend(override.aes=list(size=3)))


            figname <- paste0(volcdir, "/", compn, "_volc.pdf"); message(figname)
            pdf(figname)
            print(fig)
            bouh <- dev.off()
        }
    }


    balanced.col.grad2 <- function(mat.val) {
        min.max <- c(min(mat.val), max(mat.val))
        abs.min <- abs(min.max[1])
        max.val <- max(abs(min.max))

        rel.min <- (abs.min / max.val)
        rel.max <- (min.max[2]) / max.val

        low.col <- rgb(1 - rel.min, 1 - rel.min, 1)
        high.col <- rgb(1, 1 - rel.max, 1 - rel.max)
        

        return(c(low.col, high.col))
    }


    hm_gene_vst_expression <- function(vst.obj, gene.sel.gff, spl_label_df, fig_base, relative.values = FALSE, spl_ref = NULL, mat.gene.id = NULL, height.factor = 1, width.factor = 1) {
        ## raw VST expression of histones
        if (TRUE) {
            #   # recover histone expression data
            if (TRUE) {
                if (class(vst.obj)[1] %in% c("DESeq2", "DESeqTransform")) {
                    svsd <- vst.obj[mcols(vst.obj)$Geneid %in% gene.sel.gff$gene_id, ]
                    fgi <- factor(mcols(svsd)$Geneid, levels=gene.sel.gff$gene_id)
                    svsd <- svsd[order(fgi), ]
                    svsd <- tibble(data.frame(assay(svsd))) %>% mutate(gene_id = mcols(svsd)$Geneid) 
                } else if (!is.null(mat.gene.id)) {
                    svsd <- vst.obj[mat.gene.id %in% gene.sel.gff$gene_id, ]
                    fgi <- factor(mat.gene.id[mat.gene.id %in% gene.sel.gff$gene_id], levels = gene.sel.gff$gene_id)
                    svsd <- svsd[order(fgi), ]
                    svsd <- log10(svsd)
                    svsd <- tibble(data.frame(svsd)) %>% mutate(gene_id = fgi[order(fgi)]) 
                }

            }

            #   # build data.frame for figure
            if (TRUE) {
                fdata <- svsd %>% reshape2::melt(id.vars = "gene_id", value.name = "expression", variable.name = "sample") %>% tibble()
                fdata <- left_join(
                    fdata,
                    data.frame(mcols(gene.sel.gff))[c('gene_id', 'gene_name')],
                    by = "gene_id"
                )
                fdata$gene_name <- factor(fdata$gene_name, levels = unique(rev(gene.sel.gff$gene_name)))
                fdata <- left_join(
                    fdata,
                    spl_label_df,
                    by = "sample"
                )
                fdata$spl_label <- factor(fdata$spl_label, levels = unique(fdata$spl_label))
            }
            
            #   # build heatmap of expression
            if (TRUE) {
                fig <- ggplot(
                    fdata,
                    aes(
                        x = spl_label,
                        y = gene_name,
                        fill = expression
                    )
                ) +
                    geom_tile() +
                    
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))

                if (relative.values) {
                    col.grad2.cust <- balanced.col.grad2(fdata$expression)
                    fig <- fig + scale_fill_gradient2(low = col.grad2.cust[1], mid = "white", high = col.grad2.cust[2])
                    # fig <- fig + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
                } else {
                    fig <- fig + scale_fill_viridis_c()
                }

                # figname <- paste0(figdir, "/histone_gene_expression_hm.pdf"); message(figname)
                figname <- paste0(fig_base, "_hm.pdf"); message(figname)
                pdf(figname, height = height.factor * 2 * 7, width = width.factor * 7)
                print(fig)
                bouh <- dev.off()
            }
        }

        ## log2 fold change according to 0h of replication
        if (!is.null(spl_ref)) {
            #   # compute histone expression relative to reference sample
            if (TRUE) {
                lsvsd <- svsd
                for (colo in colnames(svsd)[grepl('^D', colnames(svsd))]) {
                    # lsvsd[[colo]] <- log2(lsvsd[[colo]] / svsd[['D356T1']])
                    lsvsd[[colo]] <- log2(lsvsd[[colo]] / svsd[[spl_ref]])
                }
            }

            #   # build data.frame for figure
            if (TRUE) {
                fdata <- lsvsd %>% reshape2::melt(id.vars = "gene_id", value.name = "log2FCto0h", variable.name = "sample") %>% tibble()
                fdata <- left_join(
                    fdata,
                    data.frame(mcols(gene.sel.gff))[c('gene_id', 'gene_name')],
                    by = "gene_id"
                )
                fdata$gene_name <- factor(fdata$gene_name, levels = unique(rev(gene.sel.gff$gene_name)))
                fdata <- left_join(
                    fdata,
                    spl_label_df,
                    by = "sample"
                )
                fdata$spl_label <- factor(fdata$spl_label, levels = unique(fdata$spl_label))
            }

            #   # build heatmap of expression
            if (TRUE) {
                col.grad2.cust <- balanced.col.grad2(fdata$log2FCto0h)
                fig <- ggplot(
                    fdata,
                    aes(
                        x = spl_label,
                        y = gene_name,
                        fill = log2FCto0h
                    )
                ) +
                    geom_tile() +
                    # scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
                    scale_fill_gradient2(low = col.grad2.cust[1], mid = "white", high = col.grad2.cust[2]) +
                    theme(axis.text.x = element_text(angle = 45, hjust=1))

                figname <- paste0(fig_base, "_relativeTo_", spl_label_df$spl_label[spl_label_df$sample == spl_ref], "_hm.pdf"); message(figname)
                pdf(figname, height = height.factor * 2 * 7, width = width.factor * 7)
                print(fig)
                bouh <- dev.off()
            }
        }
    }


    gid2name <- function(id.vec, g.tab = NULL, g.gr = NULL, gene.id.col = NULL, gene.name.col = NULL) {
        if (typeof(id.vec) != typeof(list())) {
            id.vec <- list(as.vector(id.vec))
        }

        if (is.null(g.tab)) {
            if (is.null(g.gr)) {
                message('!!! Specify a table or a GRanges object with indicated gene names and corresponding ids')
                return(NULL)
            } else {
                g.tab <- as.data.frame(mcols(g.gr))
            }
        }

        if (!is.null(gene.id.col)) {
            g.tab['gene_id'] = gtab[gene.id.col]
        }

        if (!is.null(gene.name.col)) {
            g.tab['gene_name'] = gtab[gene.name.col]
        }

        out <- lapply(id.vec, function(xxx) {
            left_join(
                tibble('gene_id' = xxx),
                g.tab,
                by = 'gene_id'
            )[['gene_name']]
        })

        return(out)
    }


    gname2id <- function(name.vec, g.tab, gene.id.col = NULL, gene.name.col = NULL) {
        if (typeof(name.vec) != typeof(list())) {
            name.vec <- list(as.vector(name.vec))
        }

        if (!is.null(gene.id.col)) {
            g.tab['gene_id'] = gtab[gene.id.col]
        }

        if (!is.null(gene.name.col)) {
            g.tab['gene_name'] = gtab[gene.name.col]
        }

        out <- lapply(name.vec, function(xxx) {
            left_join(
                tibble('gene_name' = xxx),
                g.tab,
                by = 'gene_name'
            )[['gene_id']]
        })

        return(out)
    }

}


#####################################################
#### Set up the working environment
#####################################################
if (TRUE)
{
    ## environment variables
    figdir <- "./Images/ASF1_chaperone"
    resdir <- "./Results/ASF1_chaperone"
    resdir_deseq <- paste0(resdir, "/DESeq2_results"); dir.create(resdir_deseq, showWarnings = F, recursive = F)
    read_count_cds_dir <- file.path(resdir, 'transcriptome_signal_in_SL_and_HDE/featureCounts/on_coding_sequence_of_replicative_histone_genes')
    read_count_dsSL_dir <- file.path(resdir, 'transcriptome_signal_in_SL_and_HDE/featureCounts/on_SL_downstream_region')
    read_count_exon_dir <- file.path(resdir, 'transcriptome_signal_in_SL_and_HDE/featureCounts/on_distal_exon')


    ## load list of all genes in EnsEMBL dataset
    gene.gff <- rtracklayer::import.gff('./Data/GRCh38.104.sorted.gtf')
    gene.gff <- gene.gff[mcols(gene.gff)$type == "gene"]

    ## get EnsEMBL annotation of coding sequenes for replicative histones genes
    cds_rhc_p <- './Data/ASF1_chaperone/annotations/replicative_histone_CDS.gtf'
    cdsgr <- rtracklayer::import.gff(cds_rhc_p)

    ## load the table of Spike-in metadata.
    spi.tab <- tibble(read.table('./Data/ERCC92/ERCC_Controls_Analysis.txt', h=T, sep='\t'))
    as.data.frame(head(spi.tab))


    ## load list of histone genes
    histgr <- rtracklayer::import.gff('./Data/histones.gff')
    hist_name_oldnew <- tibble(read.table('./Data/histone_geneNames_OldNew.tsv'))
    colnames(hist_name_oldnew) <- c('new', 'old')
    histgr_old <- histgr
    histgr_old$gene_name <- left_join(
        data.frame('new' = histgr$gene_name),
        hist_name_oldnew,
        by = 'new'
    )$old


    ## Document duplicated replicative histones
    {
        dup.l <- list(
            c('HIST2H4A', 'HIST2H4B'),
            c('HIST2H3C', 'HIST2H3A'),
            c('HIST2H2AA3', 'HIST2H2AA4')
        )

        dup.ll <- lapply(dup.l, function(xxx) {xxx[1]})
        names(dup.ll) <- lapply(dup.l, function(xxx) {xxx[2]})

        dupi.l <- lapply(dup.l, function(xxx) {
            left_join(
                data.frame('gene_name'=xxx),
                as.data.frame(mcols(histgr_old))[c('gene_id', 'gene_name')],
                by='gene_name'
            )$gene_id
        })

        dupi.ll <- lapply(dupi.l, function(xxx) {xxx[1]})
        names(dupi.ll) <- lapply(dupi.l, function(xxx) {xxx[2]})

    }

    ## load list of genes involved in replicative histone biogenesis
    replHistPathgr <- rtracklayer::import.gff('./Data/path_replicative_histone_synthesis.gtf')


    ## load published expression data
    enc.expr.tab.path <- './Data/ASF1_chaperone/replHist_geneExpression_in_HeLa-S3_Brooks2015RNA.tsv'
    eet <- read.csv(enc.expr.tab.path, h=T, sep='\t', comment.char='#')


    #### read counts in genes (only signal in exons)
    rcdir <- "./Data/ASF1_chaperone/RNA-Seq/read_counts"
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
        "replicate" = c(rep(1:2, 8), rep(1:2, each=2)),
        "experiment" = factor(rep(c('Sphase', 'total', 'nascent'), c(12, 4, 4)), levels=c('Sphase', 'total', 'nascent'))
    ))
    desdf$exp_comb <- paste(desdf$experiment, desdf$Sphase, sep='_')


    ## compute path to read count files
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

    ## Load the clusters by Alberto
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

    ## set of gene id for replicative histone gene in group 2 of clustering
    rhgr2id <- albcl.t$gene_id[albcl.t$cluster == 2]
    rhgr2id <- rhgr2id[rhgr2id %in% histgr_old$gene_id[grepl('^HIST', histgr_old$gene_name)]]


    ## load the Spike-in count table
    spi.cnt <- tibble(read.table(file.path(resdir, 'counts.idxstats.ERCC92.tsv'), h=T, sep='\t'))
    spi.cnt <- spi.cnt[grepl('(^D356T|ERCC)', colnames(spi.cnt))]
    spi.cnt <- spi.cnt[grepl('^ERCC', spi.cnt$ERCC),]


    ## color palette for histone genes
    color_list <- list(
        'replicative' = '#A500CD', # '#CE00FF', # '#8d7493',
        'replicative_dark' = '#67007F',
        'replacement' = '#15D900', #'#11B300' # '#7bae6c'
        'replacement_dark' = '#3D7F26'
    )

    ## dot size to highligh histone genes
    hh_size = 5
    hh_alpha = 0.80

}


#####################################################
#### Analyze the spike-in levels
#####################################################
if (TRUE)
{
    ## Normalize by sequencing depth
    {
        ## Seq depth on genome
        tot_seqdep_gen <- colSums(mcnt[2:ncol(mcnt)])
        if (SAVE_DATA_LOG) {
            save_data_f(named_vec_to_df(tot_seqdep_gen, id_name='sample', value_name='seq_depth'), resdir, 'sequencing_depth_per_sample', table=TRUE)
        }

        ## Normalization factors based on sequencing depth
        tot_seqdep_gen_nf <- min(tot_seqdep_gen) / tot_seqdep_gen
        if (SAVE_DATA_LOG) {
            save_data_f(named_vec_to_df(tot_seqdep_gen_nf, id_name='sample', value_name='seq_depth'), resdir, 'sequencing_depth_norm_factor', table=TRUE)
        }

        ## Seq depth on spike-in
        tot_seqdep_spi <- colSums(spi.cnt[2:ncol(spi.cnt)])
        tot_seqdep_spi <- tot_seqdep_spi[names(tot_seqdep_gen)]

        ## Total seq depth
        tot_seqdep <- tot_seqdep_gen + tot_seqdep_spi

        ## Get RPKM
        spi.rpkm <- apply(spi.cnt[2:ncol(spi.cnt)], 1, function(xxx) {
            (xxx / tot_seqdep_spi) * 1000
        })
        spi.rpkm <- t(spi.rpkm)

        spi.rpkm <- tibble(as.data.frame(spi.rpkm))
        # mhead(spi.rpkm)

        spi.rpkm$ERCC.ID <- spi.cnt$ERCC
        # mhead(spi.rpkm)
    }


    ## Compare expected concentration and read counts
    {
        ## build the data frame for the figure
        {
            gdata <- pivot_longer(data=spi.rpkm, cols = desdf$sample, names_to = 'sample', values_to = 'rpkm')

            gdata$rpkm <- gdata$rpkm + 0.1

            gdata$mix1 <- left_join(
                gdata['ERCC.ID'],
                spi.tab[c('ERCC.ID', 'concentration.in.Mix.1..attomoles.ul.')],
                by='ERCC.ID'
            )[['concentration.in.Mix.1..attomoles.ul.']]

            gdata$mix2 <- left_join(
                gdata['ERCC.ID'],
                spi.tab[c('ERCC.ID', 'concentration.in.Mix.2..attomoles.ul.')],
                by='ERCC.ID'
            )[['concentration.in.Mix.2..attomoles.ul.']]
        }


        ## build figure for every sample and every mix
        if (TRUE) {
            figdir_spic <- file.path(figdir, 'spikes_concentrations'); dir.create(figdir_spic, showWarnings = F, recursive = F)

            for (spln in desdf$sample) {
                for (mixn in c('mix1', 'mix2')) {
                    gdata$concentration = gdata[[mixn]]

                    ## build the figure
                    fig <- ggplot(
                        gdata[gdata$sample == spln,],
                        aes(
                            x = log10(concentration),
                            y = log10(rpkm)
                        )
                    ) +
                        geom_point() +
                        theme_txt

                    if (PRINT_PLOTS)
                    {
                        figname <- file.path(figdir_spic, paste0('spike_concentration_spl', spln, '_vs_', mixn, '.pdf')); message(figname)
                        pdf(figname)
                        print(fig)
                        bouh <- dev.off()
                    }
                }
            }
        }


        ## determine which mix used for each sample
        {
            ## threshold for minimum detection (and meaningful results)
            min_log10rpkm_thres <- -2

            ## log transform
            gdata$rpkm_log10 <- log10(gdata$rpkm)

            ## evaluate each sample for every mix
            {
                corres.l <- list()
                for (spln in desdf$sample) {
                    corres.l[[spln]] <- list()
                    for (mixn in c('mix1', 'mix2')) {
                        gdata$concentration <- gdata[[mixn]]
                        gdata$concen_log10 <- log10(gdata$concentration)
                        
                        ## get correlation RPKM vs Concentration
                        log_sel <- gdata$sample == spln & gdata$rpkm_log10 >= -2
                        corres <- cor(gdata$concen_log10[log_sel], gdata$rpkm_log10[log_sel])
                        corres.l[[spln]][[mixn]] <- corres
                    }
                }
                corres.l <- reshape2::melt(corres.l)
                colnames(corres.l) <- c('correlation', 'mix', 'sample')
                corres.l <- pivot_wider(data=corres.l, id_cols = 'sample', names_from = 'mix', values_from = 'correlation')

                corres.l$mix_good <- apply(corres.l[c('mix1', 'mix2')], 1, function(xxx) {
                    ifelse(xxx[1] > xxx[2], 'mix1', 'mix2')
                })
                corres.l$mix_good <- factor(corres.l$mix_good, levels=c('mix1', 'mix2'))
            }

            ## save the table sample / mix
            {
                if (SAVE_DATA_LOG) {
                    save_data_f(corres.l[c('sample', 'mix_good')], resdir, 'spike_in_mix_in_samples', table=TRUE)
                }
            }
        }
    }


    ## Get concentration of the Spike-in in the samples
    {
        ## compute spike-in fraction
        spi.fract <- tot_seqdep_spi / tot_seqdep

        ## get genome fraction (complementary of spike-in fraction)
        gen.fract <- 1 - spi.fract

        ## plot the spike-ing fraction
        if (TRUE) {
            fig <- ggplot(
                data.frame(
                    'fraction' = as.vector(spi.fract),
                    'sample' = factor(names(spi.fract), levels=names(spi.fract))
                ),
                aes(
                    x = sample,
                    y = fraction
                )
            ) +
                geom_col() + 
                theme_txt_xangle
            
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'spike_fraction_samples.pdf'); message(figname)
                pdf(figname)
                print(fig)
                bouh <- dev.off()
            }
        }
    }


    ## Compute the normalization factors based on spike-in
    {
        min_spike_fract <- min(spi.fract)
        spi.fract.min.rel <- min_spike_fract / spi.fract

        if (SAVE_DATA_LOG) {
            save_data_f(
                data.frame('sample' = names(spi.fract.min.rel), 'norm_factor' = as.vector(spi.fract.min.rel)),
                resdir,
                'spike_norm_factor', table=TRUE
            )
        }
    }


    ## Compute the normalization factors based on combined seq depth and spike-in
    {
        spi.seqdep.nf <- spi.fract.min.rel * tot_seqdep_gen_nf

        if (SAVE_DATA_LOG) {
            save_data_f(
                data.frame('sample' = names(spi.seqdep.nf), 'norm_factor' = as.vector(spi.seqdep.nf)),
                resdir,
                'spike_and_seqDepth_norm_factor', table=TRUE
            )
        }
    }
}


#####################################################
#### Analyze gene expression in the RNA-seq experiments
#####################################################
if (TRUE) 
{
    ## evaluate the quantitativity per gene
    {
        ## spike-normalized counts
        nmcnt <- apply(mcnt[desdf$sample], 1, function(xxx) {
            xxx * spi.seqdep.nf[desdf$sample]
        })
        nmcnt <- t(nmcnt)

        ## Repartition of genes by read counts
        nmrep <- apply(
            nmcnt,
            2,
            function(xxx) {
                repart_fc(xxx, keep.original.order=TRUE)[['repart']]
            }
        )
        nmrep <- as.data.frame(nmrep)
        nmrep$gene_id <- mcnt$Geneid

        ## complement table of spike-normalized counts
        nmcnt <- as.data.frame(nmcnt)
        nmcnt$gene_id <- mcnt$Geneid

        ## reshape and associate normalized count with repartition
        gdata <- pivot_longer(
            data = nmcnt,
            cols = desdf$sample,
            names_to = 'sample',
            values_to = 'norm_count'
        )
        gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))

        gdata$repart <- pivot_longer(
            data = nmrep,
            cols = desdf$sample,
            names_to = 'sample',
            values_to = 'repart'
        )$repart

        ## build figure
        if (TRUE)
        {
            fig <- ggplot(
                gdata[gdata$norm_count != 0,],
                aes(
                    x = log10(norm_count),
                    y = repart,
                    col = sample
                )
            ) +
                geom_point(size = 0.1) +
                geom_abline(slope = -0.025, intercept = 0.81)

            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'spike_normed_counts_repart_by_sample.png'); message(figname)
                png(figname, width = 480 * 2, height = 480 * 2)
                print(fig)
                bouh <- dev.off()
            }
        }


        ## select the valid genes
        {
            ## get the list of selections
            sgi.l <- list()
            for (expn in levels(desdf$experiment)) {
                log_sel <- gdata$norm_count > 1 & (gdata$repart - 0.81) / log10(gdata$norm_count) > -0.025 & gdata$sample %in% desdf$sample[desdf$experiment == expn]
                sel_gene_id <- unique(gdata$gene_id[log_sel])
                #
                sgi.l[[expn]] = sel_gene_id
                remove(sel_gene_id)
            }

            ## get the intersection of all the lists
            sgi.l[['common']] <- sgi.l[['total']]
            for (expn in levels(desdf$experiment)) {
                sgi.l$common <- sgi.l$common[sgi.l$common %in% sgi.l[[expn]]]
            }
            # lapply(sgi.l, length)

            ## record the selections
            if (SAVE_DATA_LOG) {
                save_data_f(sgi.l, resdir, 'spike_normed_based_gene_selection_for_DiffExprAnal')
            }
        }
    }

    

    ## Evaluate pseudo-counts and logFoldChange in 4sU- and total RNA-seq
    {
        dds.l <- list()
        dds.res.l <- list()
        for (expn in c("total", "nascent")) {
            sel_spl_vec <- desdf$sample[desdf$experiment == expn]

            ## create EdgeR object
            xxx = as.matrix(mcnt[sel_spl_vec])
            rownames(xxx) <- mcnt$Geneid
            dds <- DGEList(counts=xxx, group=desdf$condition[desdf$experiment == expn])
            # dds$samples

            ## filter genes
            # keep <- filterByExpr(dds)
            keep <- rownames(xxx) %in% sgi.l[[expn]]
            dds <- dds[keep, , keep.lib.sizes = FALSE]
            
            ## apply normalization
            dds$samples$norm.factors <- spi.seqdep.nf[sel_spl_vec]

            ## estimate common dispersion and tagwise dispersion (tag ~ gene, exon, ... whatever feature, the lines of the count matrix)
            dds <- estimateDisp(dds)
            dds.l[[expn]] <- dds

            ## Make the differential expression test
            dds.res <- exactTest(dds)

            ## distribution of p-values
            if (PRINT_PLOTS)
            {
                fig <- repart_plot(-log10(dds.res$table$PValue))

                figname <- file.path(figdir, 'test.png'); message(figname)
                png(figname)
                print(fig)
                bouh <- dev.off()
            }

            ## record the results in list
            dds.res$table$gene_id <- rownames(dds.res$table)
            dds.res$table$gene_name <- left_join(
                dds.res$table,
                tibble(as.data.frame(mcols(gene.gff)))[c('gene_id', 'gene_name')],
                by = 'gene_id'
            )$gene_name
            # str(dds.res$table$gene_name)

            dds.res.l[[expn]] = dds.res
        }

        ## save results
        if (SAVE_DATA_LOG) {
            save_data_f(dds.res.l, resdir, 'spike_normed_based_FC_Pval_and_pseudoCounts')
        }
    }


    ## Compare FoldChanges total vs nascent
    {
        ## collect the logFC infos
        {
            common_gene_ids <- dds.res.l$total$table$gene_id[dds.res.l$total$table$gene_id %in% dds.res.l$nascent$table$gene_id]
            # length(common_gene_ids)

            gdata <- dds.res.l$total$table[dds.res.l$total$table$gene_id %in% common_gene_ids, c('logFC', 'gene_id', 'gene_name')]
            gdata$logFC.total <- gdata$logFC
            gdata$logFC <- NULL
            gdata$logFC.nascent <- left_join(
                gdata['gene_id'],
                dds.res.l$nascent$table[c('gene_id', 'logFC')],
                by = 'gene_id'
            )$logFC
            # head(gdata)
        }

        ## build the figure
        {
            fig <- ggplot(
                gdata,
                aes(
                    x = logFC.total,
                    y = logFC.nascent
                )
            ) + 
                geom_point(
                    pch = 21,
                    fill = '#CCCCCC', # '#CCCCCC',
                    col = '#AAAAAA',
                    alpha = 1,
                    size = hh_size / 2
                ) + 
                geom_abline(slope=1, intercept=0, col='#AA0000') +
                geom_vline(xintercept = 0, col = '#000000') +
                geom_hline(yintercept = 0, col = '#000000') +
                geom_point(
                    data = gdata[gdata$gene_id %in% histgr_old$gene_id[!grepl('^HIST', histgr_old$gene_name)],],
                    pch=21,
                    fill = color_list$replacement,
                    col = color_list$replacement_dark,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                geom_point(
                    data = gdata[gdata$gene_id %in% rhgr2id,],
                    pch=21,
                    fill = color_list$replicative,
                    col = color_list$replicative_dark,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                theme_txt
        }

        ## print the figure
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'logFC_EdgeR_spikeNormed_total_vs_nascent_2d.pdf'); message(figname)
            pdf(figname)
            print(fig)
            bouh <- dev.off()
        }
    }


    ## Identify the replicative histones with positive logFC in total RNA-seq
    {
        dds.res.tt = dds.res.l[['total']]$table
        sel.dds.res.tt <- dds.res.tt[dds.res.tt$logFC > 0 & dds.res.tt$gene_id %in% histgr_old$gene_id[grepl('^HIST', histgr_old$gene_name)],]
    }


    ## Get log2(CPM) expression values in S-phase synchronized experiment
    {
        expn <- "Sphase"
        sel_spl_vec <- desdf$sample[desdf$experiment == expn]

        ## create EdgeR object
        xxx = as.matrix(mcnt[sel_spl_vec])
        rownames(xxx) <- mcnt$Geneid
        groups <- paste0(desdf$condition[desdf$experiment == expn], ':', desdf$Sphase[desdf$experiment == expn])
        dds <- DGEList(counts=xxx, group=factor(groups, levels=unique(groups)))
        # dds$samples

        ## filter genes
        # keep <- filterByExpr(dds)
        keep <- rownames(xxx) %in% sgi.l[[expn]]
        dds <- dds[keep, , keep.lib.sizes = FALSE]
        
        ## apply normalization
        dds$samples$norm.factors <- spi.seqdep.nf[sel_spl_vec]

        ## estimate common dispersion and tagwise dispersion (tag ~ gene, exon, ... whatever feature, the lines of the count matrix)
        dds <- estimateDisp(dds)
        dds.l[[expn]] <- dds

        logcpm <- cpm(dds, log=TRUE)
        # mhead(logcpm)

        ## compute diff expression with the control "siCtrl:0"
        for (grpi in 2:length(levels(dds$samples$group))) {
            dds.res <- exactTest(dds, pair=c(1, grpi))
            dds.res.l[[levels(dds$samples$group)[grpi]]] <- dds.res
        }

    }

    ## Identify down-regulated replicative histones
    {
        ## subset to replicative histone genes
        rhc.logcpm <- logcpm[rownames(logcpm) %in% rhgr2id,]
        # mhead(rhc.logcpm)

        ## Plot expression through the samples
        {
            ## reshape the data
            {
                gdata <- as.data.frame(rhc.logcpm)
                gdata$gene_id <- rownames(gdata)

                gdata <- pivot_longer(data = gdata, cols = desdf$sample[desdf$experiment == 'Sphase'], names_to = 'sample', values_to = 'logcpm')
                gdata$sample <- factor(gdata$sample, levels=unique(gdata$sample))
                gdata$gene_name <- left_join(
                    gdata['gene_id'],
                    as.data.frame(histgr_old)[c('gene_id', 'gene_name')],
                    by = 'gene_id'
                )$gene_name
            }


            ## build the plot
            {
                fig <- ggplot(
                    gdata,
                    aes(
                        x = sample,
                        y = logcpm,
                        col = gene_name
                    )
                ) +
                    geom_jitter(size=1, width=0.25) +
                    theme_txt_xangle
            }

            ## print the plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'spike_normed_log2CPM_replHistones_SphaseSamples.pdf'); message(figname)
                pdf(figname, width = 7 * 3)
                print(fig)
                bouh <- dev.off()
            }


            ## compute logCPM relative to control sample (siCtrl_0h)
            {
                rhc.logcpm.relctrl <- rhc.logcpm
                # mhead(rhc.logcpm)

                ## get mean logCPM in control condition
                mean_ctrl <- apply(
                    rhc.logcpm[,desdf$sample[desdf$experiment == 'Sphase' & desdf$condition == 'siCtrl' & !is.na(desdf$Sphase == 0) & desdf$Sphase == 0]],
                    1,
                    mean
                )

                ## log foldchange with this mean control logCPM
                {
                    rhc.logcpm.relctrl <- apply(
                        rhc.logcpm.relctrl,
                        2,
                        function(xxx) {
                            xxx - mean_ctrl
                        }
                    )

                    # mhead(rhc.logcpm.relctrl)
                }
            }

            ## reshape the data
            {
                hdata <- as.data.frame(rhc.logcpm.relctrl)
                hdata$gene_id <- rownames(hdata)

                hdata <- pivot_longer(data = hdata, cols = desdf$sample[desdf$experiment == 'Sphase'], names_to = 'sample', values_to = 'rel.logcpm')
                hdata$sample <- factor(hdata$sample, levels=unique(hdata$sample))
                hdata$gene_name <- left_join(
                    hdata['gene_id'],
                    as.data.frame(histgr_old)[c('gene_id', 'gene_name')],
                    by = 'gene_id'
                )$gene_name
            }

            ## build the plot
            {
                fig <- ggplot(
                    hdata,
                    aes(
                        x = sample,
                        y = rel.logcpm,
                        col = gene_name
                    )
                ) +
                    geom_jitter(size=1, width=0.25) +
                    theme_txt_xangle
            }

            ## print the plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'spike_normed_ctrlRelative_log2CPM_replHistones_SphaseSamples.pdf'); message(figname)
                pdf(figname, width = 7 * 3)
                print(fig)
                bouh <- dev.off()
            }


            ## build as a heatmap
            {
                col.pal <- rev(brewer.pal(n = 11, name = 'RdYlBu'))

                fig <- ggplot(
                    hdata %>% mutate(sample = factor(sample, levels=rev(levels(sample)))),
                    aes(
                        x = gene_name,
                        y = sample,
                        fill = rel.logcpm
                    )
                ) +
                    geom_tile() +
                    scale_fill_gradientn(colours=col.pal) +
                    theme_txt_xangle
                    # scale_fill_gradient2(low = "#3A3A98", mid = "#FFFFE9", high = '#832424', midpoint = 0, na.value = 'white') + 
            }

            ## print the heatmap
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'spike_normed_ctrlRelative_log2CPM_replHistones_SphaseSamples_hm.pdf'); message(figname)
                pdf(figname, width = 7 * 3)
                print(fig)
                bouh <- dev.off()
            }

            ## list of outliers among replicative histone genes
            {
                rhc.outl <- c('HIST1H2BA', 'HIST1H2BB', 'HIST1H2BM', 'HIST1H3C', 'HIST2H4A', 'HIST3H2BB')
            }

            ## rebuild the plot
            {
                fig <- ggplot(
                    hdata[!(hdata$gene_name %in% rhc.outl),],
                    aes(
                        x = sample,
                        y = rel.logcpm,
                        col = gene_name
                    )
                ) +
                    geom_jitter(size=1, width=0.25) +
                    theme_txt_xangle
            }

            ## print the plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'spike_normed_ctrlRelative_log2CPM_replHistones_ConsensusOnly_SphaseSamples.pdf'); message(figname)
                pdf(figname, width = 7 * 3)
                print(fig)
                bouh <- dev.off()
            }
        }
        #
    }


    ## Stratify genes by their length
    {
        ## repartition of by gene length for genes expression in 4sU-RNA-seq and total RNA-seq
        {
            gl.rpt <- repart_fc(
                width(gene.gff[gene.gff$gene_id %in% common_gene_ids]),
                keep.original.order = TRUE
            )

            gl.rpt$gene_id <- gene.gff$gene_id[gene.gff$gene_id %in% common_gene_ids]
            gl.rpt$gene_name <- gene.gff$gene_name[gene.gff$gene_id %in% common_gene_ids]
            # tibble(gl.rpt)
        }

        ## define classes of genes according to their gene length
        {
            gl.bord.log <- c(0, 2.5, 3, 3.5, 4, 5, 6)
            gl.rpt$lgtc <- 0
            for (bii in 1:length(gl.bord.log)) {
                gl.rpt$lgtc[log10(gl.rpt$values) > gl.bord.log[bii]] <- bii
            }
            gl.rpt$lgtc <- factor(as.character(gl.rpt$lgtc))
            # tibble(gl.rpt)
            # levels(gl.rpt$lgtc)
        }
        
        ## build repartition plot
        {
            # col.pal.class <- rev(brewer.pal(n = 9, name = 'Greys'))[1:7]
            col.pal.class <- hsv(h=1, s=1, v=(6:0) / 6)
            names(col.pal.class) <- levels(gl.rpt$lgtc)

            fig <- ggplot(
                gl.rpt,
                aes(
                    x = log10(values),
                    y = repart,
                    col = lgtc
                )
            ) + 
                geom_point() + 
                geom_point(
                    data = gl.rpt[gl.rpt$gene_id %in% rhgr2id,],
                    col = color_list$replicative,
                    size = hh_size
                ) + 
                scale_color_manual(values=col.pal.class) +
                theme_txt
        }

        ## print the plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'gene_length_repartition_exprIn4sUaTotal.pdf'); message(figname)
            pdf(figname)
            print(fig)
            bouh <- dev.off()
        }

        ## rebuild the 4sU vs total logFC plot with gene length classes
        {
            ## collect the logFC infos
            {
                gdata <- dds.res.l$total$table[dds.res.l$total$table$gene_id %in% common_gene_ids, c('logFC', 'gene_id', 'gene_name')]
                gdata$logFC.total <- gdata$logFC
                gdata$logFC <- NULL
                gdata$logFC.nascent <- left_join(
                    gdata['gene_id'],
                    dds.res.l$nascent$table[c('gene_id', 'logFC')],
                    by = 'gene_id'
                )$logFC
                # head(gdata)

                gdata$lgtc <-left_join(
                    gdata['gene_id'],
                    gl.rpt[c('lgtc', 'gene_id')],
                    by = 'gene_id'
                )$lgtc
                # any(is.na(gdata$lgtc))
            }

            ## build the figure
            {
                # col.pal.class <- rev(brewer.pal(n = 9, name = 'Greys'))[1:7]
                col.pal.class <- hsv(h=1, s=1, v=(6:0) / 6)
                names(col.pal.class) <- levels(gl.rpt$lgtc)

                fig <- ggplot(
                    gdata,
                    aes(
                        x = logFC.total,
                        y = logFC.nascent,
                        col = lgtc
                    )
                ) + 
                    geom_point(size=1) + #, col='#CCCCCC') +
                    geom_abline(slope=1, intercept=0, col='#AA0000') +
                    geom_vline(xintercept = 0, col = '#000000') +
                    geom_hline(yintercept = 0, col = '#000000') +
                    geom_point(
                        data = gdata[gdata$gene_id %in% histgr_old$gene_id[!grepl('^HIST', histgr_old$gene_name)],],
                        col = color_list$replacement,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    geom_point(
                        data = gdata[gdata$gene_id %in% rhgr2id,],
                        col = color_list$replicative,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    scale_color_manual(values=col.pal.class) +
                    theme_txt
            }

            ## print the figure
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'logFC_EdgeR_spikeNormed_total_vs_nascent_wGeneLengthClasses_2d.pdf'); message(figname)
                pdf(figname)
                print(fig)
                bouh <- dev.off()
            }

            ## build the figure with classes in facets
            {
                # col.pal.class <- rev(brewer.pal(n = 9, name = 'Greys'))[1:7]
                col.pal.class <- hsv(h=1, s=1, v=(6:0) / 6)
                names(col.pal.class) <- levels(gl.rpt$lgtc)

                fig <- ggplot(
                    gdata,
                    aes(
                        x = logFC.total,
                        y = logFC.nascent,
                        col = lgtc
                    )
                ) + 
                    facet_wrap(~lgtc) +
                    # geom_point(size=1) + #, col='#CCCCCC') +
                    geom_density2d() +
                    geom_abline(slope=1, intercept=0, col='#AA0000') +
                    geom_vline(xintercept = 0, col = '#000000') +
                    geom_hline(yintercept = 0, col = '#000000') +
                    geom_point(
                        data = gdata[gdata$gene_id %in% histgr_old$gene_id[!grepl('^HIST', histgr_old$gene_name)],],
                        col = color_list$replacement,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    geom_point(
                        data = gdata[gdata$gene_id %in% rhgr2id,],
                        col = color_list$replicative,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    scale_color_manual(values=col.pal.class) +
                    theme_txt
            }

            ## print the figure with facets
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'logFC_EdgeR_spikeNormed_total_vs_nascent_wGeneLengthClassesInFacets_2d.pdf'); message(figname)
                pdf(figname, width = 7*2, height = 7*2)
                print(fig)
                bouh <- dev.off()
            }
        }


        ## build scatter plot gene length vs 4sU-RNA-seq logFC
        {
            ## collect the infos
            {
                gdata <- dds.res.l$nascent$table[dds.res.l$total$table$gene_id %in% common_gene_ids, c('logFC', 'gene_id', 'gene_name')]
                gdata$gene_length <- left_join(
                    gdata,
                    gl.rpt[c('values', 'gene_id')],
                    by = 'gene_id'
                )$values

                gdata$lgtc <- left_join(
                    gdata,
                    gl.rpt[c('lgtc', 'gene_id')],
                    by = 'gene_id'
                )$lgtc

                # any(is.na(gdata$gene_length))
                # str(gdata)
            }


            ## build the figure
            {
                col.pal.class <- hsv(h=1, s=1, v=(6:0) / 6)
                names(col.pal.class) <- levels(gl.rpt$lgtc)

                fig <- ggplot(
                    gdata,
                    aes(
                        x = log10(gene_length),
                        y = logFC,
                        col = lgtc
                    )
                ) + 
                    geom_point(size=1, col='#CCCCCC') +
                    geom_point(
                        data = gdata[gdata$gene_id %in% histgr_old$gene_id[!grepl('^HIST', histgr_old$gene_name)],],
                        pch=21,
                        fill = color_list$replacement,
                        col = color_list$replacement_dark,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    geom_point(
                        data = gdata[gdata$gene_id %in% rhgr2id,],
                        pch=21,
                        fill = color_list$replicative,
                        col = color_list$replicative_dark,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    scale_color_manual(values=col.pal.class) +
                    theme_txt
            }

            ## print the figure
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'logFC_EdgeR_spikeNormed_nascent_vs_GeneLength.pdf'); message(figname)
                pdf(figname)
                print(fig)
                bouh <- dev.off()
            }

        }
    }
}



#####################################################
#### Focus on the replicative histones CDSs
#####################################################
if (TRUE)
{
    ## Load Read Counts in CDS of replicative histone genes
    {
        ## Compile the counts of reads per replicative histone
        rtab <- list()
        seq_depth_list <- list()
        for (filename in list.files(read_count_cds_dir)) {
            # if (grepl('D356T[0-9]*_wholeGene.rv$', filename)) {
            if (grepl('replicative_histone_CDS_D356T[0-9]*.tsv$', filename)) {
                message(filename)
                intab <- read.table(paste0(read_count_cds_dir, '/', filename), h = T)
                seq_depth_list[[filename]] <- sum(intab[[7]])
                intab <- intab[intab$Geneid %in% cdsgr$gene_id,]
                rtab[[filename]] <- intab
            }
        }

        #### Compile the total mRNA read counts for the replication histone genes
        if (TRUE) {
            ## reduce the count tables in one unique table
            rcdata <- tibble(rtab[[1]])
            rcdata <- rcdata[1:ncol(rcdata) - 1]
            rcdata['repl_histone'] <- left_join(
                rcdata['Geneid'] %>% mutate(gene_id = Geneid),
                as.data.frame(mcols(gene.gff))[c('gene_id', 'gene_name')],
                by = 'gene_id'
            )[['gene_name']]


            for (name in names(rtab)) {
                spl_name <- strsplit(strsplit(name, '\\.')[[1]][1], '_')[[1]][4]
                rcdata[[spl_name]] <- rtab[[name]][[7]]
            }
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

        ## select the replicative histone genes from group 2
        rcdata <- rcdata[rcdata$Geneid %in% rhgr2id,]

        ## normalize by sequencing depth and spike-in
        {
            for (spln in desdf$sample) {
                rcdata[[spln]] <- rcdata[[spln]] * spi.seqdep.nf[spln]
            }
        }

    }


    ## Evaluate pseudo-counts and logFoldChange in 4sU- and total RNA-seq
    {
        cds.dds.l <- list()
        cds.dds.res.l <- list()
        for (expn in c("total", "nascent")) {
            # expn <- 'nascent'
            sel_spl_vec <- desdf$sample[desdf$experiment == expn]

            ## create EdgeR object
            xxx = as.matrix(rcdata[sel_spl_vec])
            rownames(xxx) <- rcdata$Geneid
            cds.dds <- DGEList(counts=xxx, group=desdf$condition[desdf$experiment == expn])
            # cds.dds$samples

            ## filter genes
            # keep <- filterByExpr(cds.dds)
            keep <- rownames(xxx) %in% rcdata$Geneid
            cds.dds <- cds.dds[keep, , keep.lib.sizes = FALSE]
            
            ## apply normalization
            cds.dds$samples$norm.factors <- spi.seqdep.nf[sel_spl_vec]

            ## estimate common dispersion and tagwise dispersion (tag ~ gene, exon, ... whatever feature, the lines of the count matrix)
            cds.dds <- estimateDisp(cds.dds)
            cds.dds.l[[expn]] <- cds.dds

            ## Make the differential expression test
            cds.dds.res <- exactTest(cds.dds)

            ## record the results in list
            cds.dds.res$table$gene_id <- rownames(cds.dds.res$table)
            cds.dds.res$table$gene_name <- left_join(
                cds.dds.res$table,
                tibble(as.data.frame(mcols(gene.gff)))[c('gene_id', 'gene_name')],
                by = 'gene_id'
            )$gene_name
            # str(cds.dds.res$table$gene_name)

            cds.dds.res.l[[expn]] = cds.dds.res
        }

        ## save results
        if (SAVE_DATA_LOG) {
            save_data_f(cds.dds.res.l, resdir, 'spike_normed_based_FC_Pval_and_pseudoCounts_ReplHistGeneCDS')
        }
    }


    ## Compare logFC and logCPM in whole gene vs CDS annotations
    {
        for (expn in c('total', 'nascent')) {
            # expn <- 'nascent'
            message(expn)

            ## recover the data
            {
                cds.dds.res <- cds.dds.res.l[[expn]]
                
                gdata <- cds.dds.res$table[c('gene_id', 'gene_name', 'logFC', 'logCPM')]
                gdata$logFC.CDS <- gdata$logFC
                gdata$logCPM.CDS <- gdata$logCPM

                gdata$logFC <- NULL
                gdata$logCPM <- NULL

                gdata$logFC.tot <- left_join(
                    gdata,
                    dds.res.l[[expn]]$table,
                    by = 'gene_id'
                )$logFC

                gdata$logCPM.tot <- left_join(
                    gdata,
                    dds.res.l[[expn]]$table,
                    by = 'gene_id'
                )$logCPM

                gdata <- gdata[!is.na(gdata$logFC.tot),]
            }

            ## focus on replicative histone in group 2
            {
                gdata <- gdata[gdata$gene_id %in% rhgr2id,]
            }

            ## compute a regression line
            {
                out_gene_name <- c('H4C15', 'H3C14', 'H2BC6', 'H3C12')

                reg.line <- lm(logFC.tot ~ logFC.CDS, gdata[!(gdata$gene_name %in% out_gene_name),])
                print(reg.line)
            }


            ## build figure
            {
                fig <- ggplot(
                    gdata,
                    aes(
                        x = logFC.CDS,
                        y = logFC.tot
                    )
                ) +
                    geom_abline(slope = reg.line$coefficients[2], intercept = reg.line$coefficients[1], col='#CC0000') +
                    geom_point(
                        pch=21,
                        fill = color_list$replicative,
                        col = color_list$replicative_dark,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    theme_txt
            }


            ## print the figure 
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, paste0('spike_normed_wholeGeneVsCDS_replHistone_logFC_', expn, '.pdf')); message(figname)
                pdf(figname)
                print(fig)
                bouh <- dev.off()
            }
        }
    }
}


#####################################################
#### Balance CDS vs Downstream region of SL
#####################################################
if (TRUE)
{
    ## Load Read Counts in downstream region from SL of replicative histone genes
    {
        ## Compile the counts of reads per replicative histone
        dssl.rtab <- list()
        dssl.seq_depth_list <- list()
        for (filename in list.files(read_count_dsSL_dir)) {
            if (grepl('^D356T[0-9]*_replicative_histone_downstream_of_SL_region\\.rv$', filename)) {
                message(filename)
                intab <- read.table(paste0(read_count_dsSL_dir, '/', filename), h = T)
                intab$Geneid <- left_join(
                    intab %>% mutate(gene_name = Geneid),
                    as.data.frame(mcols(gene.gff)),
                    by='gene_name'
                )$gene_id
                dssl.seq_depth_list[[filename]] <- sum(intab[[7]])
                intab <- intab[intab$Geneid %in% cdsgr$gene_id,]
                dssl.rtab[[filename]] <- intab
            }
        }



        #### Compile the total mRNA read counts for the replication histone genes
        {
            ## reduce the count tables in one unique table
            dssl.rcdata <- tibble(dssl.rtab[[1]])
            dssl.rcdata <- dssl.rcdata[1:ncol(dssl.rcdata) - 1]
            dssl.rcdata['repl_histone'] <- left_join(
                dssl.rcdata['Geneid'] %>% mutate(gene_id = Geneid),
                as.data.frame(mcols(gene.gff))[c('gene_id', 'gene_name')],
                by = 'gene_id'
            )[['gene_name']]


            for (name in names(dssl.rtab)) {
                spl_name <- strsplit(strsplit(name, '\\.')[[1]][1], '_')[[1]][1]
                dssl.rcdata[[spl_name]] <- dssl.rtab[[name]][[7]]
            }
        }


        ## combine the duplicated replicative histone genes
        for (dgi in 1:length(dupi.l)) {
            # dgi <- 1
            dgp <- dupi.l[[dgi]]
            
            sdssl.rcdata <- dssl.rcdata[dssl.rcdata$Geneid %in% dgp,]
            # print(stab)

            sel_col <- !grepl('^(Geneid|Chr|Start|End|Strand|Length|repl_histone)', colnames(sdssl.rcdata))
            sdssl.rcdata[1, sel_col] <- sdssl.rcdata[1, sel_col] + sdssl.rcdata[2, sel_col]

            dssl.rcdata[dssl.rcdata$Geneid %in% dgp[1], sel_col] <- sdssl.rcdata[1, sel_col]
            dssl.rcdata <- dssl.rcdata[!(dssl.rcdata$Geneid %in% dgp[2]),]
        }

        ## select the replicative histone genes from group 2
        dssl.rcdata <- dssl.rcdata[dssl.rcdata$Geneid %in% rhgr2id,]


        ## normalize by sequencing depth and spike-in
        {
            for (spln in desdf$sample) {
                dssl.rcdata[[spln]] <- dssl.rcdata[[spln]] * spi.seqdep.nf[spln]
            }
        }

        # head(dssl.rcdata)
    }


    ## Balance CDS vs SL-downstream region
    {
        ## recover normalized read counts for various plots
        {
            ## get data for SL-downstream region
            gdata <- dssl.rcdata[grepl('^(D356T|Geneid)', colnames(dssl.rcdata))]
            gdata <- pivot_longer(
                gdata,
                cols=colnames(gdata)[grepl('^D356T', colnames(gdata))],
                values_to = 'count_dsSL',
                names_to = 'sample'
            )

            ## get data for CDS region
            pdata <- rcdata[grepl('^(D356T|Geneid)', colnames(rcdata))]
            pdata <- pivot_longer(
                pdata,
                cols=colnames(pdata)[grepl('^D356T', colnames(pdata))],
                values_to = 'count_CDS',
                names_to = 'sample'
            )

            ## compile the two
            gdata$count_CDS <- left_join(
                gdata %>% mutate(id=paste0(Geneid, '_', sample)),
                (pdata %>% mutate(id=paste0(Geneid, '_', sample)))[c('id', 'count_CDS')],
                by = 'id'
            )$count_CDS

            gdata$condition <- left_join(
                gdata,
                desdf,
                by='sample'
            )$condition

            # gdata$experiment <- left_join(
            #     gdata,
            #     desdf,
            #     by='sample'
            # )$experiment

            gdata$exp_comb <- left_join(
                gdata,
                desdf,
                by='sample'
            )$exp_comb
        }

        ## save counts in dsSL and CDS regions
        {
            outdf <- gdata
            outdf$gene_name <- left_join(
                as.data.frame(outdf %>% mutate(gene_id=Geneid)),
                as.data.frame(mcols(gene.gff)[c('gene_id', 'gene_name')]),
                by = 'gene_id'
            )$gene_name

            outdf$ratio_dsSL_over_CDS <- outdf$count_dsSL / outdf$count_CDS
            
            outdf <- outdf[c("Geneid", "gene_name", "sample", "condition", "exp_comb", "count_dsSL", "count_CDS", "ratio_dsSL_over_CDS")]
            save_data_f(outdf, resdir, 'countsAndRatio_dsSL_over_CDS', table = TRUE)

        }


        ## some heatmap plots specific to Total RNA-seq signal
        {
            ## build heatmap for raw signal in 4sU-RNA-seq
            {
                sgdata <- gdata[gdata$sample %in% desdf$sample[desdf$experiment == 'total'],]
                sgdata <- pivot_longer(sgdata, cols=c("count_dsSL", "count_CDS"), names_to = 'region', values_to = 'signal')
                sgdata$hm_col <- paste0(sgdata$region, '_', sgdata$condition, '_', sgdata$sample)
                sgdata$hm_col <- factor(sgdata$hm_col, levels=c("count_CDS_siCtrl_D356T13", "count_CDS_siCtrl_D356T14", "count_CDS_siAsf1_D356T15", "count_CDS_siAsf1_D356T16", "count_dsSL_siCtrl_D356T13", "count_dsSL_siCtrl_D356T14", "count_dsSL_siAsf1_D356T15", "count_dsSL_siAsf1_D356T16"))

                fig <- ggplot(
                    sgdata,
                    aes(
                        x = hm_col,
                        y = Geneid,
                        fill = log10(signal)
                    )
                ) +
                    geom_tile() +
                    scale_fill_viridis_c() +
                    theme_txt_xangle
            }


            ## print heatmap plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_total_hm.pdf'); message(figname)
                pdf(figname, width = 7)
                print(fig)
                bouh <- dev.off()
            }


            ## build heatmap for siCtrl_CDS-normalized signal in 4sU-RNA-seq
            {
                custfc <- function(sgdata) {
                    rsgdata = sgdata[c("Geneid", "hm_col", "signal")]

                    ## reshape with a column per 'region'
                    tmat <- pivot_wider(
                        sgdata[c("Geneid", "hm_col", "signal")],
                        names_from = 'hm_col',
                        values_from = 'signal'
                    )
                    
                    ## normalize the counts by those of CDS in siCtrl condition
                    {
                        tmat$count_CDS_siAsf1_D356T15 <- tmat$count_CDS_siAsf1_D356T15 / tmat$count_CDS_siCtrl_D356T13
                        tmat$count_dsSL_siCtrl_D356T13 <- tmat$count_dsSL_siCtrl_D356T13 / tmat$count_CDS_siCtrl_D356T13
                        tmat$count_dsSL_siAsf1_D356T15 <- tmat$count_dsSL_siAsf1_D356T15 / tmat$count_CDS_siCtrl_D356T13

                        tmat$count_CDS_siAsf1_D356T16 <- tmat$count_CDS_siAsf1_D356T16 / tmat$count_CDS_siCtrl_D356T14
                        tmat$count_dsSL_siCtrl_D356T14 <- tmat$count_dsSL_siCtrl_D356T14 / tmat$count_CDS_siCtrl_D356T14
                        tmat$count_dsSL_siAsf1_D356T16 <- tmat$count_dsSL_siAsf1_D356T16 / tmat$count_CDS_siCtrl_D356T14

                        tmat$count_CDS_siCtrl_D356T13 <- tmat$count_CDS_siCtrl_D356T13 / tmat$count_CDS_siCtrl_D356T13
                        tmat$count_CDS_siCtrl_D356T14 <- tmat$count_CDS_siCtrl_D356T14 / tmat$count_CDS_siCtrl_D356T14
                    }

                    ## get the row indexes in the original table
                    {
                        imat <- pivot_wider(
                            sgdata[c("Geneid", "hm_col")] %>% add_column(idx=1:nrow(sgdata)),
                            names_from = 'hm_col',
                            values_from = 'idx'
                        )
                        imat <- as.vector(as.matrix(imat[2:ncol(imat)]))
                    }

                    ## replace the values
                    {
                        rsgdata$signal[imat] <- as.vector(as.matrix(tmat[2:ncol(tmat)]))
                    }

                    return(rsgdata)
                }
                rsgdata <- custfc(sgdata)

                fig <- ggplot(
                    rsgdata,
                    aes(
                        x = hm_col,
                        y = Geneid,
                        fill = log10(signal)
                    )
                ) +
                    geom_tile() +
                    scale_fill_gradient2(low='blue', high='red') +
                    theme_txt_xangle
            }


            ## print heatmap plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_normBySiCtrlCDS_total_hm.pdf'); message(figname)
                pdf(figname, width = 7)
                print(fig)
                bouh <- dev.off()
            }


            ## build heatmap for siCtrl-normalized signal in 4sU-RNA-seq, CDS and 3' of SL separately
            {
                custfc <- function(sgdata) {
                    rsgdata = sgdata[c("Geneid", "hm_col", "signal")]

                    ## reshape with a column per 'region'
                    tmat <- pivot_wider(
                        sgdata[c("Geneid", "hm_col", "signal")],
                        names_from = 'hm_col',
                        values_from = 'signal'
                    )
                    
                    ## normalize the counts by those of CDS in siCtrl condition
                    {
                        tmat$count_CDS_siAsf1_D356T15 <- tmat$count_CDS_siAsf1_D356T15 / tmat$count_CDS_siCtrl_D356T13
                        tmat$count_dsSL_siAsf1_D356T15 <- tmat$count_dsSL_siAsf1_D356T15 / tmat$count_dsSL_siCtrl_D356T13

                        tmat$count_CDS_siAsf1_D356T16 <- tmat$count_CDS_siAsf1_D356T16 / tmat$count_CDS_siCtrl_D356T14
                        tmat$count_dsSL_siAsf1_D356T16 <- tmat$count_dsSL_siAsf1_D356T16 / tmat$count_dsSL_siCtrl_D356T14
                    }

                    ## get the row indexes in the original table
                    {
                        imat <- pivot_wider(
                            sgdata[c("Geneid", "hm_col")] %>% add_column(idx=1:nrow(sgdata)),
                            names_from = 'hm_col',
                            values_from = 'idx'
                        )
                        imat <- as.vector(as.matrix(imat[2:ncol(imat)]))
                    }

                    ## replace the values
                    {
                        rsgdata$signal[imat] <- as.vector(as.matrix(tmat[2:ncol(tmat)]))
                    }

                    return(rsgdata)
                }
                rsgdata <- custfc(sgdata)
                rsgdata <- rsgdata[rsgdata$hm_col %in% c("count_CDS_siAsf1_D356T15", "count_dsSL_siAsf1_D356T15", "count_CDS_siAsf1_D356T16", "count_dsSL_siAsf1_D356T16"),]

                fig <- ggplot(
                    rsgdata,
                    aes(
                        x = hm_col,
                        y = Geneid,
                        fill = log10(signal)
                    )
                ) +
                    geom_tile() +
                    scale_fill_gradient2(low='blue', high='red') +
                    theme_txt_xangle
            }


            ## print heatmap plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_normBySiCtrl_CDSandDsSL_total_hm.pdf'); message(figname)
                pdf(figname, width = 7)
                print(fig)
                bouh <- dev.off()
            }
        }


        ## some heatmap plots specific 4sU-RNA-seq signal
        {
            ## build heatmap for raw signal in 4sU-RNA-seq
            {
                sgdata <- gdata[gdata$sample %in% desdf$sample[desdf$experiment == 'nascent'],]
                sgdata <- pivot_longer(sgdata, cols=c("count_dsSL", "count_CDS"), names_to = 'region', values_to = 'signal')
                sgdata$hm_col <- paste0(sgdata$region, '_', sgdata$condition, '_', sgdata$sample)
                sgdata$hm_col <- factor(sgdata$hm_col, levels=c("count_CDS_siCtrl_D356T17", "count_CDS_siCtrl_D356T19", "count_CDS_siAsf1_D356T18", "count_CDS_siAsf1_D356T20", "count_dsSL_siCtrl_D356T17", "count_dsSL_siCtrl_D356T19", "count_dsSL_siAsf1_D356T18", "count_dsSL_siAsf1_D356T20"))

                fig <- ggplot(
                    sgdata,
                    aes(
                        x = hm_col,
                        y = Geneid,
                        fill = log10(signal)
                    )
                ) +
                    geom_tile() +
                    scale_fill_viridis_c() +
                    theme_txt_xangle
            }


            ## print heatmap plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_nascent_hm.pdf'); message(figname)
                pdf(figname, width = 7)
                print(fig)
                bouh <- dev.off()
            }


            ## build heatmap for siCtrl_CDS-normalized signal in 4sU-RNA-seq
            {
                custfc <- function(sgdata) {
                    rsgdata = sgdata[c("Geneid", "hm_col", "signal")]

                    ## reshape with a column per 'region'
                    tmat <- pivot_wider(
                        sgdata[c("Geneid", "hm_col", "signal")],
                        names_from = 'hm_col',
                        values_from = 'signal'
                    )
                    
                    ## normalize the counts by those of CDS in siCtrl condition
                    {
                        tmat$count_CDS_siAsf1_D356T18 <- tmat$count_CDS_siAsf1_D356T18 / tmat$count_CDS_siCtrl_D356T17
                        tmat$count_dsSL_siCtrl_D356T17 <- tmat$count_dsSL_siCtrl_D356T17 / tmat$count_CDS_siCtrl_D356T17
                        tmat$count_dsSL_siAsf1_D356T18 <- tmat$count_dsSL_siAsf1_D356T18 / tmat$count_CDS_siCtrl_D356T17

                        tmat$count_CDS_siAsf1_D356T20 <- tmat$count_CDS_siAsf1_D356T20 / tmat$count_CDS_siCtrl_D356T19
                        tmat$count_dsSL_siCtrl_D356T19 <- tmat$count_dsSL_siCtrl_D356T19 / tmat$count_CDS_siCtrl_D356T19
                        tmat$count_dsSL_siAsf1_D356T20 <- tmat$count_dsSL_siAsf1_D356T20 / tmat$count_CDS_siCtrl_D356T19

                        tmat$count_CDS_siCtrl_D356T17 <- tmat$count_CDS_siCtrl_D356T17 / tmat$count_CDS_siCtrl_D356T17
                        tmat$count_CDS_siCtrl_D356T19 <- tmat$count_CDS_siCtrl_D356T19 / tmat$count_CDS_siCtrl_D356T19
                    }

                    ## get the row indexes in the original table
                    {
                        imat <- pivot_wider(
                            sgdata[c("Geneid", "hm_col")] %>% add_column(idx=1:nrow(sgdata)),
                            names_from = 'hm_col',
                            values_from = 'idx'
                        )
                        imat <- as.vector(as.matrix(imat[2:ncol(imat)]))
                    }

                    ## replace the values
                    {
                        rsgdata$signal[imat] <- as.vector(as.matrix(tmat[2:ncol(tmat)]))
                    }

                    return(rsgdata)
                }
                rsgdata <- custfc(sgdata)

                fig <- ggplot(
                    rsgdata,
                    aes(
                        x = hm_col,
                        y = Geneid,
                        fill = log10(signal)
                    )
                ) +
                    geom_tile() +
                    scale_fill_gradient2(low='blue', high='red') +
                    theme_txt_xangle
            }


            ## print heatmap plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_normBySiCtrlCDS_nascent_hm.pdf'); message(figname)
                pdf(figname, width = 7)
                print(fig)
                bouh <- dev.off()
            }


            ## build heatmap for siCtrl-normalized signal in 4sU-RNA-seq, CDS and 3' of SL separately
            {
                custfc <- function(sgdata) {
                    rsgdata = sgdata[c("Geneid", "hm_col", "signal")]

                    ## reshape with a column per 'region'
                    tmat <- pivot_wider(
                        sgdata[c("Geneid", "hm_col", "signal")],
                        names_from = 'hm_col',
                        values_from = 'signal'
                    )
                    
                    ## normalize the counts by those of CDS in siCtrl condition
                    {
                        tmat$count_CDS_siAsf1_D356T18 <- tmat$count_CDS_siAsf1_D356T18 / tmat$count_CDS_siCtrl_D356T17
                        tmat$count_dsSL_siAsf1_D356T18 <- tmat$count_dsSL_siAsf1_D356T18 / tmat$count_dsSL_siCtrl_D356T17

                        tmat$count_CDS_siAsf1_D356T20 <- tmat$count_CDS_siAsf1_D356T20 / tmat$count_CDS_siCtrl_D356T19
                        tmat$count_dsSL_siAsf1_D356T20 <- tmat$count_dsSL_siAsf1_D356T20 / tmat$count_dsSL_siCtrl_D356T19
                    }

                    ## get the row indexes in the original table
                    {
                        imat <- pivot_wider(
                            sgdata[c("Geneid", "hm_col")] %>% add_column(idx=1:nrow(sgdata)),
                            names_from = 'hm_col',
                            values_from = 'idx'
                        )
                        imat <- as.vector(as.matrix(imat[2:ncol(imat)]))
                    }

                    ## replace the values
                    {
                        rsgdata$signal[imat] <- as.vector(as.matrix(tmat[2:ncol(tmat)]))
                    }

                    return(rsgdata)
                }
                rsgdata <- custfc(sgdata)
                rsgdata <- rsgdata[rsgdata$hm_col %in% c("count_CDS_siAsf1_D356T18", "count_dsSL_siAsf1_D356T18", "count_CDS_siAsf1_D356T20", "count_dsSL_siAsf1_D356T20"),]

                fig <- ggplot(
                    rsgdata,
                    aes(
                        x = hm_col,
                        y = Geneid,
                        fill = log10(signal)
                    )
                ) +
                    geom_tile() +
                    scale_fill_gradient2(low='blue', high='red') +
                    theme_txt_xangle
            }


            ## print heatmap plot
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_normBySiCtrl_CDSandDsSL_nascent_hm.pdf'); message(figname)
                pdf(figname, width = 7)
                print(fig)
                bouh <- dev.off()
            }
        }


        ## build scatter plot
        {
            fig <- ggplot(
                gdata,
                aes(
                    x = log10(count_CDS + 0.1),
                    y = log10(count_dsSL + 0.1),
                    col = condition
                )
            ) +
                # facet_wrap(~experiment, nrow=1) +
                facet_wrap(~exp_comb, nrow=1) +
                geom_point() +
                scale_color_manual(values=c('siCtrl' = '#000000', 'siAsf1' = '#CC0000')) +
                theme_txt
        }


        ## print scatter plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## build 2D density plot
        {
            fig <- ggplot(
                gdata,
                aes(
                    x = log10(count_CDS + 0.1),
                    y = log10(count_dsSL + 0.1),
                    col = condition
                )
            ) +
                facet_wrap(~exp_comb, nrow=1) +
                geom_density2d() +
                scale_color_manual(values=c('siCtrl' = '#000000', 'siAsf1' = '#CC0000')) +
                theme_txt
        }


        ## print scatter plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_dens2d.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## recover data for segment plot
        {
            ## vector of combined samples (no distinction between replicates)
            spl_comb <- paste(desdf$experiment, desdf$condition, desdf$Sphase, sep='_')

            ## get data for SL-downstream region
            agdata <- dssl.rcdata[grepl('^(D356T|Geneid)', colnames(dssl.rcdata))]
            for (spl_combn in spl_comb) {
                spls <- desdf$sample[spl_comb == spl_combn]
                agdata[[spl_combn]] <- rowMeans(agdata[spls])
            }
            agdata <- agdata[!(colnames(agdata) %in% desdf$sample)]
            
            agdata <- pivot_longer(
                agdata, 
                cols=colnames(agdata)[!grepl('^Geneid', colnames(agdata))],
                values_to = 'mean_count_dsSL',
                names_to = 'spl_comb'
            )

            agdata$condition <- factor(do.call(rbind, strsplit(agdata$spl_comb, split='_'))[, 2], levels=levels(desdf$condition))
            agdata$spl_meta <- apply(do.call(rbind, strsplit(agdata$spl_comb, split='_'))[, c(1, 3)], 1, function(xxx) {paste0(xxx[1], '_', xxx[2])})

            agdata <- pivot_wider(
                agdata,
                id_cols=c('Geneid', 'spl_meta'),
                names_from = 'condition',
                values_from = 'mean_count_dsSL'
            )
            agdata <- mutate(agdata, siCtrl_dsSL = siCtrl, siCtrl = NULL, siAsf1_dsSL = siAsf1, siAsf1 = NULL)
            # as.data.frame(agdata)


            ## get data for CDS region
            apdata <- rcdata[grepl('^(D356T|Geneid)', colnames(rcdata))]
            for (spl_combn in spl_comb) {
                spls <- desdf$sample[spl_comb == spl_combn]
                apdata[[spl_combn]] <- rowMeans(apdata[spls])
            }
            apdata <- apdata[!(colnames(apdata) %in% desdf$sample)]

            apdata <- pivot_longer(
                apdata,
                cols=colnames(apdata)[!grepl('^Geneid', colnames(apdata))],
                values_to = 'mean_count_CDS',
                names_to = 'spl_comb'
            )

            apdata$condition <- factor(do.call(rbind, strsplit(apdata$spl_comb, split='_'))[, 2], levels=levels(desdf$condition))
            apdata$spl_meta <- apply(do.call(rbind, strsplit(apdata$spl_comb, split='_'))[, c(1, 3)], 1, function(xxx) {paste0(xxx[1], '_', xxx[2])})

            apdata <- pivot_wider(
                apdata,
                id_cols=c('Geneid', 'spl_meta'),
                names_from = 'condition',
                values_from = 'mean_count_CDS'
            )
            apdata <- mutate(apdata, siCtrl_CDS = siCtrl, siCtrl = NULL, siAsf1_CDS = siAsf1, siAsf1 = NULL)
            # as.data.frame(apdata)


            ## compile the two measures
            agdata[c('siCtrl_CDS', 'siAsf1_CDS')] <- left_join(
                agdata %>% mutate(id = paste0(Geneid, '_', spl_meta)),
                apdata[c('Geneid', 'spl_meta', 'siCtrl_CDS', 'siAsf1_CDS')] %>% mutate(id = paste0(Geneid, '_', spl_meta)),
                by = 'id'
            )[c('siCtrl_CDS', 'siAsf1_CDS')]

            agdata$experiment <- do.call(rbind, strsplit(agdata$spl_meta, split='_'))[, 1]
        }


        ## build segment plot
        {
            fig <- ggplot(
                agdata,
                aes(
                    x = log10(siCtrl_CDS + 0.1),
                    xend = log10(siAsf1_CDS + 0.1),
                    y = log10(siCtrl_dsSL + 0.1),
                    yend = log10(siAsf1_dsSL + 0.1)
                )
            ) +
                # facet_wrap(~experiment, nrow=1) + 
                facet_wrap(~spl_meta, nrow=1) + 
                geom_segment(alpha=0.8) +
                theme_txt
        }


        ## print segment plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_segments.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## compute siCtrl to siAsf1 ratio displacement
        {
            agdata$log_CDS_dplcmt <- log10(agdata$siAsf1_CDS + 0.1) - log10(agdata$siCtrl_CDS + 0.1)
            agdata$log_dsSL_dplcmt <- log10(agdata$siAsf1_dsSL + 0.1) - log10(agdata$siCtrl_dsSL + 0.1)
        }


        ## Compute statistics on CDS and on SL
        {
            for (spl_meta in unique(agdata$spl_meta)) {
                # spl_meta <- unique(agdata$spl_meta)[1]
                message(paste0('>>> ', spl_meta))

                sagdata <- agdata[agdata$spl_meta == spl_meta,]

                ## in CDS
                xxx <- sagdata$log_CDS_dplcmt
                rrr <- xxx - mean(xxx)
                message('In CDS')
                t.res <- t.test(xxx, rrr, alternative='two.sided')
                print(t.res)

                ## in dsSL
                xxx <- sagdata$log_dsSL_dplcmt
                rrr <- xxx - mean(xxx)
                message('In dsSL')
                t.res <- t.test(xxx, rrr, alternative='two.sided')
                print(t.res)
            }
        }


        ## regression line for nascent experiment
        {
            out_dplcmnt <- c('ENSG00000273703', 'ENSG00000203812', 'ENSG00000272196', 'ENSG00000270882')
            sagdata <- agdata[grepl('^nascent', agdata$spl_meta) & !(agdata$Geneid %in% out_dplcmnt),]
            reg.line <- lm(
                log_dsSL_dplcmt ~ log_CDS_dplcmt,
                sagdata
            )
            print(reg.line)
        }

        ## test the slope of the regression line
        {
            ## regression
            ragdata <- sagdata
            ragdata$log_dsSL_dplcmt <- ragdata$log_dsSL_dplcmt / reg.line$coefficients[2]
            reg.line.ref <- lm(
                log_dsSL_dplcmt ~ log_CDS_dplcmt,
                ragdata
            )

            ## comparison between slopes of regression lines
            require(lsmeans)
            #   # join the observation and the computed reference
            magdata <- rbind(sagdata, ragdata) %>% mutate(lmv = c(rep('obs', nrow(sagdata)), rep('ref', nrow(ragdata))))

            #   # compute the linear regression in respect of the type of data
            m.interaction <- lm(log_dsSL_dplcmt ~ log_CDS_dplcmt*lmv, data = magdata)
            # anova(m.interaction)

            #   # compute the significance between the two slopes
            m.lst <- lstrends(m.interaction, "lmv", var="log_CDS_dplcmt")
            print(pairs(m.lst))
        }


        ## build scatter plot of displacement
        {
            fig <- ggplot(
                agdata,
                aes(
                    x = log_CDS_dplcmt,
                    y = log_dsSL_dplcmt,
                )
            ) +
                # facet_wrap(~experiment, nrow=1) + 
                facet_wrap(~spl_meta, nrow=1) + 
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_abline(slope=1, intercept = 0, col='#888888', linewidth = 3) +
                geom_point(
                    pch=21,
                    fill = color_list$replicative,
                    col = color_list$replicative_dark,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                geom_abline(
                    data = data.frame(slope=reg.line$coefficients[2], intercept=reg.line$coefficients[1], spl_meta='nascent_NA'),
                    mapping = aes(
                        slope = slope,
                        intercept = intercept
                    ),
                    col = '#CC0000',
                    linewidth = 3
                ) +
                theme_txt
        }


        ## print scatter plot of displacement
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_dplcmt2d.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## focus on "total" data
        {
            ## build scatter plot of displacement
            {
                fig <- ggplot(
                    agdata[agdata$experiment == 'total',],
                    aes(
                        x = log_CDS_dplcmt,
                        y = log_dsSL_dplcmt,
                    )
                ) +
                    # facet_wrap(~experiment, nrow=1) + 
                    # facet_wrap(~spl_meta, nrow=1) + 
                    geom_vline(xintercept = 0) +
                    geom_hline(yintercept = 0) +
                    geom_point(
                        pch=21,
                        fill = color_list$replicative,
                        col = color_list$replicative_dark,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    coord_cartesian(xlim=c(-0.7, 1), ylim=c(-0.7, 1)) +
                    theme_txt
            }


            ## print scatter plot of displacement
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_total_dplcmt2d.pdf'); message(figname)
                    pdf(figname)
                    print(fig)
                    bouh <- dev.off()
            }
        }

        ## focus on "nascent" data
        {
            ## build scatter plot of displacement
            {
                fig <- ggplot(
                    agdata[agdata$experiment == 'nascent',],
                    aes(
                        x = log_CDS_dplcmt,
                        y = log_dsSL_dplcmt,
                    )
                ) +
                    # facet_wrap(~experiment, nrow=1) + 
                    # facet_wrap(~spl_meta, nrow=1) + 
                    geom_vline(xintercept = 0) +
                    geom_hline(yintercept = 0) +
                    geom_abline(slope=1, intercept = 0, col='#888888', linewidth = 1) +
                    geom_point(
                        pch=21,
                        fill = color_list$replicative,
                        col = color_list$replicative_dark,
                        size = hh_size,
                        alpha = hh_alpha
                    ) +
                    geom_abline(
                        data = data.frame(slope=reg.line$coefficients[2], intercept=reg.line$coefficients[1], spl_meta='nascent_NA'),
                        mapping = aes(
                            slope = slope,
                            intercept = intercept
                        ),
                        col = '#CC0000',
                        linewidth = 1
                    ) +
                    coord_cartesian(xlim=c(-0.7, 1), ylim=c(-0.7, 1)) +
                    theme_txt
            }


            ## print scatter plot of displacement
            if (PRINT_PLOTS)
            {
                figname <- file.path(figdir, 'replicative_histone_CDS_vs_SLdownstream_readCounts_nascent_dplcmt2d.pdf'); message(figname)
                    pdf(figname)
                    print(fig)
                    bouh <- dev.off()
            }
        }
    }
}


#####################################################
#### Balance CDS vs distal Exon
#####################################################
if (TRUE)
{
    ## Load Read Counts in the distal exon of some of replicative histone genes
    {
        ## Compile the counts of reads per replicative histone
        ex.rtab <- list()
        ex.seq_depth_list <- list()
        for (filename in list.files(read_count_exon_dir)) {
            if (grepl('^D356T[0-9]*_replicative_histone_exons\\.rv$', filename)) {
                message(filename)
                intab <- read.table(paste0(read_count_exon_dir, '/', filename), h = T)
                intab$Geneid <- left_join(
                    intab %>% mutate(gene_name = Geneid),
                    as.data.frame(mcols(gene.gff)),
                    by='gene_name'
                )$gene_id
                ex.seq_depth_list[[filename]] <- sum(intab[[7]])
                intab <- intab[intab$Geneid %in% cdsgr$gene_id,]
                ex.rtab[[filename]] <- intab
            }
        }

        #### Compile the total mRNA read counts for the replication histone genes
        {
            ## reduce the count tables in one unique table
            ex.rcdata <- tibble(ex.rtab[[1]])
            ex.rcdata <- ex.rcdata[1:ncol(ex.rcdata) - 1]
            ex.rcdata['repl_histone'] <- left_join(
                ex.rcdata['Geneid'] %>% mutate(gene_id = Geneid),
                as.data.frame(mcols(gene.gff))[c('gene_id', 'gene_name')],
                by = 'gene_id'
            )[['gene_name']]


            for (name in names(ex.rtab)) {
                spl_name <- strsplit(strsplit(name, '\\.')[[1]][1], '_')[[1]][1]
                ex.rcdata[[spl_name]] <- ex.rtab[[name]][[7]]
            }
        }

        ## normalize by sequencing depth and spike-in
        {
            for (spln in desdf$sample) {
                ex.rcdata[[spln]] <- ex.rcdata[[spln]] * spi.seqdep.nf[spln]
            }
        }

        # head(ex.rcdata)
    }


    ## Balance CDS vs distal exon region
    {
        ## recover data for scatter plot
        {
            ## get data for distal exon region
            gdata <- ex.rcdata[grepl('^(D356T|Geneid)', colnames(ex.rcdata))]
            gdata <- pivot_longer(
                gdata,
                cols=colnames(gdata)[grepl('^D356T', colnames(gdata))],
                values_to = 'count_exon',
                names_to = 'sample'
            )

            ## get data for CDS region
            pdata <- rcdata[grepl('^(D356T|Geneid)', colnames(rcdata))]
            pdata <- pivot_longer(
                pdata,
                cols=colnames(pdata)[grepl('^D356T', colnames(pdata))],
                values_to = 'count_CDS',
                names_to = 'sample'
            )

            ## compile the two
            gdata$count_CDS <- left_join(
                gdata %>% mutate(id=paste0(Geneid, '_', sample)),
                (pdata %>% mutate(id=paste0(Geneid, '_', sample)))[c('id', 'count_CDS')],
                by = 'id'
            )$count_CDS

            gdata$condition <- left_join(
                gdata,
                desdf,
                by='sample'
            )$condition

            # gdata$experiment <- left_join(
            #     gdata,
            #     desdf,
            #     by='sample'
            # )$experiment

            gdata$exp_comb <- left_join(
                gdata,
                desdf,
                by='sample'
            )$exp_comb
        }


        ## build scatter plot
        {
            fig <- ggplot(
                gdata,
                aes(
                    x = log10(count_CDS + 0.1),
                    y = log10(count_exon + 0.1),
                    col = condition
                )
            ) +
                # facet_wrap(~experiment, nrow=1) +
                facet_wrap(~exp_comb, nrow=1) +
                geom_point() +
                scale_color_manual(values=c('siCtrl' = '#000000', 'siAsf1' = '#CC0000')) +
                theme_txt
        }


        ## print scatter plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_distalExon_readCounts.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## build 2D density plot
        {
            fig <- ggplot(
                gdata,
                aes(
                    x = log10(count_CDS + 0.1),
                    y = log10(count_exon + 0.1),
                    col = condition
                )
            ) +
                facet_wrap(~exp_comb, nrow=1) +
                geom_density2d() +
                scale_color_manual(values=c('siCtrl' = '#000000', 'siAsf1' = '#CC0000')) +
                theme_txt
        }


        ## print scatter plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_distalExon_readCounts_dens2d.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## recover data for segment plot
        {
            ## vector of combined samples (no distinction between replicates)
            spl_comb <- paste(desdf$experiment, desdf$condition, desdf$Sphase, sep='_')

            ## get data for SL-downstream region
            agdata <- ex.rcdata[grepl('^(D356T|Geneid)', colnames(ex.rcdata))]
            for (spl_combn in spl_comb) {
                spls <- desdf$sample[spl_comb == spl_combn]
                agdata[[spl_combn]] <- rowMeans(agdata[spls])
            }
            agdata <- agdata[!(colnames(agdata) %in% desdf$sample)]
            
            agdata <- pivot_longer(
                agdata, 
                cols=colnames(agdata)[!grepl('^Geneid', colnames(agdata))],
                values_to = 'mean_count_exon',
                names_to = 'spl_comb'
            )

            agdata$condition <- factor(do.call(rbind, strsplit(agdata$spl_comb, split='_'))[, 2], levels=levels(desdf$condition))
            agdata$spl_meta <- apply(do.call(rbind, strsplit(agdata$spl_comb, split='_'))[, c(1, 3)], 1, function(xxx) {paste0(xxx[1], '_', xxx[2])})

            agdata <- pivot_wider(
                agdata,
                id_cols=c('Geneid', 'spl_meta'),
                names_from = 'condition',
                values_from = 'mean_count_exon'
            )
            agdata <- mutate(agdata, siCtrl_exon = siCtrl, siCtrl = NULL, siAsf1_exon = siAsf1, siAsf1 = NULL)
            # as.data.frame(agdata)


            ## get data for CDS region
            apdata <- rcdata[grepl('^(D356T|Geneid)', colnames(rcdata))]
            for (spl_combn in spl_comb) {
                spls <- desdf$sample[spl_comb == spl_combn]
                apdata[[spl_combn]] <- rowMeans(apdata[spls])
            }
            apdata <- apdata[!(colnames(apdata) %in% desdf$sample)]

            apdata <- pivot_longer(
                apdata,
                cols=colnames(apdata)[!grepl('^Geneid', colnames(apdata))],
                values_to = 'mean_count_CDS',
                names_to = 'spl_comb'
            )

            apdata$condition <- factor(do.call(rbind, strsplit(apdata$spl_comb, split='_'))[, 2], levels=levels(desdf$condition))
            apdata$spl_meta <- apply(do.call(rbind, strsplit(apdata$spl_comb, split='_'))[, c(1, 3)], 1, function(xxx) {paste0(xxx[1], '_', xxx[2])})

            apdata <- pivot_wider(
                apdata,
                id_cols=c('Geneid', 'spl_meta'),
                names_from = 'condition',
                values_from = 'mean_count_CDS'
            )
            apdata <- mutate(apdata, siCtrl_CDS = siCtrl, siCtrl = NULL, siAsf1_CDS = siAsf1, siAsf1 = NULL)
            # as.data.frame(apdata)


            ## compile the two measures
            agdata[c('siCtrl_CDS', 'siAsf1_CDS')] <- left_join(
                agdata %>% mutate(id = paste0(Geneid, '_', spl_meta)),
                apdata[c('Geneid', 'spl_meta', 'siCtrl_CDS', 'siAsf1_CDS')] %>% mutate(id = paste0(Geneid, '_', spl_meta)),
                by = 'id'
            )[c('siCtrl_CDS', 'siAsf1_CDS')]

            agdata$experiment <- do.call(rbind, strsplit(agdata$spl_meta, split='_'))[, 1]
        }


        ## build segment plot
        {
            fig <- ggplot(
                agdata,
                aes(
                    x = log10(siCtrl_CDS + 0.1),
                    xend = log10(siAsf1_CDS + 0.1),
                    y = log10(siCtrl_exon + 0.1),
                    yend = log10(siAsf1_exon + 0.1)
                )
            ) +
                # facet_wrap(~experiment, nrow=1) + 
                facet_wrap(~spl_meta, nrow=1) + 
                geom_segment(alpha=0.8) +
                theme_txt
        }


        ## print segment plot
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_distalExon_readCounts_segments.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }


        ## compute siCtrl to siAsf1 ratio displacement
        {
            agdata$log_CDS_dplcmt <- log10(agdata$siAsf1_CDS + 0.1) - log10(agdata$siCtrl_CDS + 0.1)
            agdata$log_exon_dplcmt <- log10(agdata$siAsf1_exon + 0.1) - log10(agdata$siCtrl_exon + 0.1)
        }


        ## build scatter plot of displacement
        {
            fig <- ggplot(
                agdata,
                aes(
                    x = log_CDS_dplcmt,
                    y = log_exon_dplcmt,
                )
            ) +
                # facet_wrap(~experiment, nrow=1) + 
                facet_wrap(~spl_meta, nrow=1) + 
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_abline(slope=1, intercept = 0, col='#888888') +
                geom_point(
                    pch=21,
                    fill = color_list$replicative,
                    col = color_list$replicative_dark,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                theme_txt
        }

        ## print scatter plot of displacement
        if (PRINT_PLOTS)
        {
            figname <- file.path(figdir, 'replicative_histone_CDS_vs_distalExon_readCounts_dplcmt2d.pdf'); message(figname)
            pdf(figname, width = 5 * 7)
            print(fig)
            bouh <- dev.off()
        }
    }
}




#####################################################
#### Compare expression level with published data
#####################################################
if (TRUE)
{
    ## combine data from published study and from our study
    {
        gdata <- eet[c('gene_name', 'ENCODE_HeLa.S3_polyA._FPKM')]

        ## retrieve replicative histone gene expression from our samples
        {
            log_rhc <- mcnt$Geneid %in% histgr_old$gene_id
            rhc.nmcnt <- nmcnt[log_rhc,]

            rhc.nmcnt$gene_name <- left_join(
                rhc.nmcnt['gene_id'],
                as.data.frame(mcols(histgr_old))[c('gene_name', 'gene_id')],
                by = 'gene_id'
            )$gene_name
        }

        ## Summarize the duplicated genes
        {
            for (iii in 1:length(dup.l)) {
                # iii <- 1
                tmp <- apply(rhc.nmcnt[rhc.nmcnt$gene_name %in% dup.l[[iii]], desdf$sample], 2, sum)
                rhc.nmcnt[rhc.nmcnt$gene_name == dup.l[[iii]][1], desdf$sample] <- tmp
                rhc.nmcnt <- rhc.nmcnt[!(rhc.nmcnt$gene_name == dup.l[[iii]][2]),]
            }
        }

        ## gather all in one table
        gdata <- left_join(
            gdata,
            rhc.nmcnt,
            by='gene_name'
        )
    }


    ## build the figures
    {
        fig.l <- list()
        for (spln in desdf$sample) {
            pdata <- gdata
            pdata$sample <- pdata[[spln]]
            #
            fig.l[[spln]] <- ggplot(
                pdata,
                aes(
                    x = log10(ENCODE_HeLa.S3_polyA._FPKM + 0.1),
                    y = log10(sample + 0.1)
                )
            ) + 
                geom_point() +
                theme_txt
        }
    }


    ## print the figures
    if (PRINT_PLOTS)
    {
        figdir_compexpr <- file.path(figdir, 'ours_vs_published_HeLa_expr'); dir.create(figdir_compexpr, showWarnings = F, recursive = F)

        for (spln in desdf$sample) {
            figname <- file.path(figdir_compexpr, paste0(spln, '_vs_published.pdf')); message(figname)
            pdf(figname)
            print(fig.l[[spln]])
            bouh <- dev.off()
        }
    }
}



#####################################################
#### Analyze remapped 'polyA' reads
#####################################################
if (TRUE)
{
    # require(Rsamtools)
    

    resdir_polya <- file.path(resdir, 'poly_adenylation') #, 'fq_reads_wPolyAremoved')

    ## load the identified polyA-reads
    {
        samt.l <- list()
        samgi.l <- list()
        for (spln in desdf$sample) {
            samf <- file.path(resdir_polya, paste0('replicative_histone_genes_', spln, '_polyAreads.sam'))
            samt <- read.table(samf, h=F, sep='\n', comment.char = "@")
            samt <- as.list(samt[[1]])

            samt <- lapply(samt, function(xxx) {
                strsplit(xxx, split='\t')[[1]]
            })

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
    }


    ## get contingencies of each genes polyA reads
    {
        allgi <- unique(unlist(samgi.l))
        samgim <- matrix(0, nrow=length(allgi), ncol=length(desdf$sample), dimnames=list(allgi, desdf$sample))
        samgi.cnt <- lapply(samgi.l, table)
        for (spln in desdf$sample) {
            for (gin in names(samgi.cnt[[spln]])) {
                samgim[gin, spln] <- samgi.cnt[[spln]][gin]
            }
        }
    }


    ## build gene tile distribution plot
    {
        gdata <- pivot_longer(data.frame(samgim) %>% mutate(gene_id=rownames(samgim)), desdf$sample, names_to = 'sample', values_to = 'count')
        gdata$count[gdata$count == 0] <- NA
        gdata$sample <- factor(gdata$sample, levels=desdf$sample)

        gdata$gene_name <- left_join(
            gdata,
            as.data.frame(mcols(histgr_old))[c('gene_id', 'gene_name')],
            by = 'gene_id'
        )$gene_name

        fig <- ggplot(
            gdata,
            aes(
                x = sample,
                y = gene_name,
                fill = count
            )
        ) +
            geom_tile() +
            scale_fill_gradient2(midpoint=0.5,) +
            theme_txt_xangle
    }


    ## print plot
    if (PRINT_PLOTS)
    {
        figname <- file.path(figdir, 'polyA', 'gene_occurence.pdf'); message(figname)
        pdf(figname, width = 2 * 7)
        print(fig)
        dev.off()
    }


    ## look for PAS 'AAUAAA' in the polyA read sequences
    {
        has.pas.nbr.l <- list()
        for (spln in names(samt.l)) {
            samt <- samt.l[[spln]]
            
            has.pas <- lapply(samt, function(xxx) {
                grepl('TTTATT', xxx[10])
            })

            has.pas.nbr <- sum(unlist(has.pas))
            has.pas.nbr.l[[spln]] <- has.pas.nbr
        }
        has.pas.nbr.l <- named_vec_to_df(unlist(has.pas.nbr.l))
        colnames(has.pas.nbr.l) <- c('sample', 'count')
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
        gdata$group <- paste0(gdata$exp_comb, '_', gdata$condition)
        gdata$group <- factor(gdata$group, levels=unique(gdata$group))

        fig <- ggplot(
            gdata,
            aes(
                x = group,
                y = count,
                group = group
            )
        ) +
            facet_wrap(~experiment, scales = "free_x") +
            geom_point() + 
            geom_segment(
                data = pivot_wider(gdata %>% mutate(repl_names = paste0('repl', replicate)), id_cols=c('group', 'experiment'), names_from = 'repl_names', values_from = 'count'),
                mapping = aes(
                    x = group,
                    y = repl1,
                    xend = group,
                    yend = repl2
                )
            ) +
            theme_txt_xangle
    }

    ## print the plot
    if (PRINT_PLOTS)
    {
        figname <- file.path(figdir, 'polyA', 'AAUAAA_occurences_in_polyAreads.pdf'); message(figname)
        pdf(figname)
        print(fig)
        bouh <- dev.off()
    }
}

####
