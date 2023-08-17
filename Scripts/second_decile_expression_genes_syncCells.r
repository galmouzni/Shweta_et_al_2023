if (TRUE) {
    require(reshape2)
    require(data.table)
    require(tidyverse)
    require(ggdendro)
    require(grid)
    require(DESeq2)
    library("vsn")
    library("pheatmap")
    library("RColorBrewer")
    library("IHW")
    require(igraph)
    require(WebGestaltR)
    library("biomaRt")
}

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
if (TRUE) {
    ## environment variables
    figdir <- "./Images/ASF1_chaperone"
    resdir <- "./Results/ASF1_chaperone"
    resdir_deseq <- paste0(resdir, "/DESeq2_results"); dir.create(resdir_deseq, showWarnings = F, recursive = F)

    ## load list of all genes in EnsEMBL dataset
    gene.gff <- rtracklayer::import.gff('./Data/GRCh38.104.sorted.gtf')
    gene.gff <- gene.gff[mcols(gene.gff)$type == "gene"]


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

    ## load list of genes involved in replicative histone biogenesis
    replHistPathgr <- rtracklayer::import.gff('./Data/path_replicative_histone_synthesis.gtf')

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
        "replicate" = c(rep(1:2, 8), rep(1:2, each=2))
    ))

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
}

#####################################################
#### on synchronized cells at 0h, 2h, 5h post G1/S arrest; in conditions +/- Asf1
#####################################################
if (TRUE) {
    spl_label_df <- data.frame(
        "sample" = c("D356T1", "D356T2", "D356T3", "D356T4", "D356T5", "D356T6", "D356T7", "D356T8", "D356T9", "D356T10", "D356T11", "D356T12"),
        "spl_label" = c("ctrl.0h.1", "ctrl.0h.2", "ctrl.2h.1", "ctrl.2h.2", "ctrl.5h.1", "ctrl.5h.2", "siASF1.0h.1", "siASF1.0h.2", "siASF1.2h.1", "siASF1.2h.2", "siASF1.5h.1", "siASF1.5h.2")
    )

    gene_id_vec <- mcnt[, 'Geneid']

    ## select samples of sync cells RNA-seq
    smcnt <- mcnt[,colnames(mcnt) %in% spl_label_df$sample]

    ## filter out genes not expressed throughout the RNA-seq
    log_expr <- rowSums(smcnt) > 10
    smcnt <- smcnt[log_expr,]
    expr_gene_id_vec <- gene_id_vec[log_expr,]

    ## compute order of genes by increasing expression in each samples
    smrk <- smcnt
    for (spl in colnames(smcnt)) {
        smrk[[spl]] <- order(order(smcnt[[spl]]))
    }

    ## compute average order
    avgrk <- rowMeans(smrk)
    
    ## get decile for each gene
    gdec_vec <- ceiling((avgrk / length(avgrk)) * 10)
    
    ## compute decile group
    dec_list <- list()
    for (deci in 1:10) {
        dec_list[[deci]] <- as.vector(unlist(expr_gene_id_vec[gdec_vec == deci,]))
    }

    ## export GFF per decile
    resdir_decile <- file.path(resdir, 'decile_expr_genes_increasing_syncCellsRNAseq'); dir.create(resdir_decile, showWarnings = F, recursive = F)
    for (deci in 1:10) {
        dec_gff <- gene.gff[gene.gff$gene_id %in% dec_list[[deci]]]
        rtracklayer::export.gff3(dec_gff, file.path(resdir_decile, paste0('decile_', deci, '.gff')))
    }
}


####
