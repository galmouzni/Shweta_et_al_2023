## Differential gene expression analysis of Schweta's project (+/- ASF1)

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

oldtheme <- theme_set(theme_bw())

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

if (TRUE) {
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }


    df_to_named_vec <- function(df, value_col='value', name_col='name') {
        ## convert a vector of named item into a dataframe
        nvec <- df[[value_col]]
        names(nvec) <- df[[name_col]]

        return(nvec)
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
}


#####################################################
#### analysis of S-phase time points, +/- Asf1
#####################################################
if (TRUE) {
    spl_label_df <- data.frame(
        "sample" = c("D356T1", "D356T2", "D356T3", "D356T4", "D356T5", "D356T6", "D356T7", "D356T8", "D356T9", "D356T10", "D356T11", "D356T12"),
        "spl_label" = c("ctrl.0h.1", "ctrl.0h.2", "ctrl.2h.1", "ctrl.2h.2", "ctrl.5h.1", "ctrl.5h.2", "siASF1.0h.1", "siASF1.0h.2", "siASF1.2h.1", "siASF1.2h.2", "siASF1.5h.1", "siASF1.5h.2")
    )

    ## DESeq2 setting
    dds <- DESeqDataSetFromMatrix(countData = mcnt[,2:13],
                                colData = desdf[1:12,],
                                design = ~ condition + Sphase + Sphase:condition)


    mcols(dds) <- DataFrame(mcols(dds), data.frame(Geneid = mcnt["Geneid"]))

    dds$condition <- factor(dds$condition, levels = c("siCtrl", "siAsf1"))

    ## independent filtering
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]


    ## DESeq processing
    dds <- DESeq(dds)

    ## transforming data (VST)
    vsd <- vst(dds, blind=FALSE)

    ## check independency between mean / standard deviation
    figname <- paste0(figdir, '/meanSdPlot_VST.pdf'); message(figname)
    pdf(figname)
    meanSdPlot(assay(vsd))
    bouh <- dev.off()



    ## distance heatmap/hclust of the samples
    sampleDists <- dist(t(assay(vsd))) # default: method = "euclidean"
    sampleDistMatrix <- as.matrix(sampleDists)

    rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$Sphase, sep=":")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors)


    ## exploring dataset structure
    plotPCA(vsd, intgroup=c("condition", "Sphase"))


    # resultsNames(dds)
    resl <- list(
        "Asf1VsCtrl" = results(dds, contrast=c("condition", "siAsf1", "siCtrl")), # main effect of siAsf1 vs siCtrl on 0h of S-phase
        "2Vs0" = results(dds, contrast=c("Sphase", "2", "0")), # main effect of S-phase at 2h vs 0h in condition of siCtrl
        "5Vs0" = results(dds, contrast=c("Sphase", "5", "0")), # main effect of S-phase at 5h vs 0h in condition of siCtrl
        "5Vs2" = results(dds, contrast=c("Sphase", "5", "2")), # main effect of S-phase at 5h vs 2h in condition of siCtrl
        "5From2Vs0" = results(dds, contrast=list(c("Sphase_5_vs_0", "Sphase_2_vs_0"))), # additional effect of S-phase at 5h to effet of 2h vs 0h in condition of siCtrl
        "A2VsC0" = results(dds, name="conditionsiAsf1.Sphase2"), # test if the effec of S-phase at 2h vs 0h is different in condition siAsf1 compared to siCtrl
        "A5VsC0" = results(dds, name="conditionsiAsf1.Sphase5"), # test if the effec of S-phase at 5h vs 0h is different in condition siAsf1 compared to siCtrl
        "A5FromA2VsC0" = results(dds, contrast=list(c("conditionsiAsf1.Sphase5", "conditionsiAsf1.Sphase2"))) # variation between 5h compared to 2h vs 0h of S-phase in the additional effect by siAsf1 vs siCtrl.
    )

    if (FALSE) {
        resIHW <- results(dds, filterFun=ihw)
        summary(resIHW)
        metadata(resIHW)$ihwResult
    }

    ##  # MA plot
    if (TRUE) {
        par(mfrow=c(2,4))
        for (compn in names(resl)) {
            plotMA(resl[[compn]], alpha = 0.01)
            title(compn)
        }
        par(mfrow=c(1,1))
    }


    ## build volcano plot for every comparison
    volcdir <- paste0(figdir, "/CtrlVsSiASF1_Sphase_comparisonVolcanoPlot")
    dir.create(volcdir, showWarnings = F, recursive = F)
    
    volcano_plots(resl, volcdir)


    ##  # look at the expression of the DEGs by comparison
    if (TRUE) {
        figname <- paste0(figdir, '/', 'DEG_expression_by_comparison_boxplot.pdf')
        pdf(figname, width = 4 * 7, height = 7)
        nres = length(resl)
        grid.newpage()
        idx <- 0
        for (degt in c("up", "down")) {
            for (compn in names(resl)) {
                dtab <- vsd[(!is.na(resl[[compn]]$padj)) & resl[[compn]]$padj <= 0.01 & list("up" = resl[[compn]]$log2FoldChange > 0, "down" = resl[[compn]]$log2FoldChange < 0)[[degt]]]
                gdata <- tibble(gather(as.data.frame(assay(dtab)), sample, expression))
                gdata <- inner_join(gdata, desdf, by = "sample")
                gdata$sample_ <- paste0(gdata$condition, '_', gdata$Sphase)
                gdata$spl_name <- paste0(gdata$sample_, '_', gdata$replicate)
                gdata$spl_name <- factor(gdata$spl_name, levels = unique(gdata$spl_name), order = T)
                
                fig <- ggplot(gdata, aes(x = Sphase, y = expression, fill = condition, group = spl_name)) + 
                    geom_boxplot() + 
                    scale_fill_manual(values = gg_color_hue(2)) + 
                    ggtitle(paste0(compn, '_', degt)) + 
                    theme(legend.position = 'bottom')

                print(fig, vp = viewport(x = 1/(nres*2) + (idx %% nres) * 1/nres, y = 1 - (1/(2*2) + (idx %/% nres) * 1/2), width = 1/nres, height = 1/2))
                idx <- idx + 1
            }
        }
        bouh <- dev.off()
    }


    ## get summary about the whole set of DEGs
    if (TRUE) {
        fdr.thresh <- 0.01
        deg.summary <- deg.summary.fc(resl=resl, expr.mtx=vsd, fdr.thresh = fdr.thresh)

        ## save the summary
        sfn <- paste0(resdir, '/ctrlVsSiASF1_Sphase_DEG_VSTexpr_dendrogram_summary.RDS')
        saveRDS(object = deg.summary, file = sfn)
    }
}



#####################################################
#### analysis of +/- Asf1 in asynchronous cells
#####################################################
if (TRUE) {
    asynch.total.log <- is.na(desdf$Sphase) & desdf$RNA == "total"

    ## DESeq2 setting
    at.dds <- DESeqDataSetFromMatrix(countData = mcnt[,which(asynch.total.log) + 1],
                                colData = desdf[asynch.total.log,],
                                design = ~ condition)


    mcols(at.dds) <- DataFrame(mcols(at.dds), data.frame(Geneid = mcnt["Geneid"]))
    mcols(at.dds)['GeneOS'] <- left_join(
        data.frame(mcols(at.dds))['Geneid'],
        data.frame(mcols(gene.gff))[c('gene_id', 'gene_name')] %>% dplyr::mutate(Geneid = gene_id),
        by = 'Geneid'
    )[['gene_name']]

    at.dds$condition <- factor(at.dds$condition, levels = c("siCtrl", "siAsf1"))

    ## independent filtering
    keep <- rowSums(counts(at.dds)) >= 10
    at.dds <- at.dds[keep,]


    ## DESeq processing
    at.dds <- DESeq(at.dds)

    ## transforming data (VST)
    at.vsd <- vst(at.dds, blind=FALSE)


    ## distance heatmap/hclust of the samples
    sampleDists <- dist(t(assay(at.vsd))) # default: method = "euclidean"
    sampleDistMatrix <- as.matrix(sampleDists)

    rownames(sampleDistMatrix) <- at.vsd$condition
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    figname <- paste0(figdir, "/Asynchronous_siASF1_effect_sample_correlation_hclust.pdf"); message(figname)
    pdf(figname)
    print(
        pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors)
    )
    bouh <- dev.off()


    ## exploring dataset structure
    figname <- paste0(figdir, "/Asynchronous_siASF1_effect_sample_PCA.pdf"); message(figname)
    pdf(figname)
    print(
        plotPCA(at.vsd, intgroup=c("condition"))
    )
    bouh <- dev.off()


    # resultsNames(dds)
    at.resl <- list(
        "Asf1VsCtrl" = results(at.dds, contrast=c("condition", "siAsf1", "siCtrl")) # main effect of siAsf1 vs siCtrl on asynchronous cells
    )


    ## build volcano plot for every comparison
    at.volcdir <- paste0(figdir, "/CtrlVsSiASF1_Async_comparisonVolcanoPlot")
    dir.create(at.volcdir, showWarnings = F, recursive = F)
    
    volcano_plots(at.resl, at.volcdir)


    ##  # look at the expression of the DEGs by comparison
    fdr.thresh <- 0.01
    if (TRUE) {
        figname <- paste0(figdir, '/', 'DEG_expression_onAsync_by_comparison_boxplot.pdf')
        pdf(figname, width = 0.25 * 7, height = 7)
        nres = length(at.resl)
        grid.newpage()
        idx <- 0
        for (degt in c("up", "down")) {
            for (compn in names(at.resl)) {
                dtab <- at.vsd[(!is.na(at.resl[[compn]]$padj)) & at.resl[[compn]]$padj <= fdr.thresh & list("up" = at.resl[[compn]]$log2FoldChange > 0, "down" = at.resl[[compn]]$log2FoldChange < 0)[[degt]]]
                gdata <- tibble(gather(as.data.frame(assay(dtab)), sample, expression))
                gdata <- inner_join(gdata, desdf, by = "sample")
                gdata$sample_ <- paste0(gdata$condition)
                gdata$spl_name <- paste0(gdata$sample_, '_', gdata$replicate)
                gdata$spl_name <- factor(gdata$spl_name, levels = unique(gdata$spl_name), order = T)
                
                fig <- ggplot(gdata, aes(x = spl_name, y = expression, fill = condition, group = spl_name)) + 
                    geom_boxplot() + 
                    scale_fill_manual(values = gg_color_hue(2)) + 
                    ggtitle(paste0(compn, '_', degt)) + 
                    theme(legend.position = 'bottom')

                print(fig, vp = viewport(x = 1/(nres*2) + (idx %% nres) * 1/nres, y = 1 - (1/(2*2) + (idx %/% nres) * 1/2), width = 1/nres, height = 1/2))
                idx <- idx + 1
            }
        }
        bouh <- dev.off()
    }


    ## get summary about the whole set of DEGs
    if (TRUE) {
        at.deg.summary <- deg.summary.fc(resl=at.resl, expr.mtx=at.vsd, fdr.thresh = fdr.thresh)

        ## save the summary
        at.sfn <- paste0(resdir, '/ctrlVsSiASF1_Async_DEG_VSTexpr_dendrogram_summary.RDS')
        saveRDS(object = at.deg.summary, file = at.sfn)
    }
}


#####################################################
#### analysis of S-phase situation compared to asynchronous (~average) situation
#####################################################
if (TRUE) {
    ctrl.total.log <- desdf$condition == "siCtrl" & desdf$RNA == "total"

    ## DESeq2 setting
    ctrl.dds <- DESeqDataSetFromMatrix(countData = mcnt[,which(ctrl.total.log) + 1],
                                colData = tibble(desdf[ctrl.total.log,]) %>% mutate(Sphase = ifelse(is.na(Sphase), 'none', as.character(Sphase))),
                                design = ~ Sphase)
    ctrl.dds@colData@listData$Sphase <- factor(ctrl.dds@colData@listData$Sphase, levels=c('none', '0', '2', '5'))

    mcols(ctrl.dds) <- DataFrame(mcols(ctrl.dds), data.frame(Geneid = mcnt["Geneid"]))

    ## independent filtering
    keep <- rowSums(counts(ctrl.dds)) >= 10
    ctrl.dds <- ctrl.dds[keep,]


    ## DESeq processing
    ctrl.dds <- DESeq(ctrl.dds)

    ## transforming data (VST)
    ctrl.vsd <- vst(ctrl.dds, blind=FALSE)


    ## distance heatmap/hclust of the samples
    sampleDists <- dist(t(assay(ctrl.vsd))) # default: method = "euclidean"
    sampleDistMatrix <- as.matrix(sampleDists)

    rownames(sampleDistMatrix) <- ctrl.vsd$Sphase
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    figname <- paste0(figdir, "/Ctrl_AsyncVsSphase_effect_sample_correlation_hclust.pdf"); message(figname)
    pdf(figname)
    print(
        pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors)
    )
    bouh <- dev.off()


    ## exploring dataset structure
    figname <- paste0(figdir, "/Ctrl_AsyncVsSphase_effect_sample_PCA.pdf"); message(figname)
    pdf(figname)
    print(
        plotPCA(ctrl.vsd, intgroup=c("Sphase"))
    )
    bouh <- dev.off()


    # resultsNames(dds)
    ctrl.resl <- list(
        "AsyncVs0" = results(ctrl.dds, contrast=c("Sphase", "0", "none")),
        "AsyncVs2" = results(ctrl.dds, contrast=c("Sphase", "2", "none")),
        "AsyncVs5" = results(ctrl.dds, contrast=c("Sphase", "5", "none"))
    )


    ##  # look at the expression of the DEGs by comparison
    fdr.thresh <- 0.01
    l2fc.thresh <- 1.5 # only used for this chunck
    if (TRUE) {
        figname <- paste0(figdir, '/', 'DEG_expression_AsyncVsSphase_by_comparison_boxplot.pdf')
        pdf(figname, width = 1 * 7, height = 7)
        nres = length(ctrl.resl)
        grid.newpage()
        idx <- 0
        for (degt in c("up", "down")) {
            for (compn in names(ctrl.resl)) {
                dtab <- ctrl.vsd[(!is.na(ctrl.resl[[compn]]$padj)) & ctrl.resl[[compn]]$padj <= fdr.thresh & abs(ctrl.resl[[compn]]$log2FoldChange) >= l2fc.thresh & list("up" = ctrl.resl[[compn]]$log2FoldChange > 0, "down" = ctrl.resl[[compn]]$log2FoldChange < 0)[[degt]]]
                gdata <- tibble(gather(as.data.frame(assay(dtab)), sample, expression))
                gdata <- inner_join(gdata, desdf, by = "sample")
                gdata$sample_ <- paste0(gdata$condition, '_', gdata$Sphase)
                gdata$spl_name <- paste0(gdata$sample_, '_', gdata$replicate)
                gdata$spl_name <- factor(gdata$spl_name, levels = unique(gdata$spl_name), order = T)
                
                fig <- ggplot(gdata, aes(x = Sphase, y = expression, fill = condition, group = spl_name)) + 
                    geom_boxplot() + 
                    scale_fill_manual(values = gg_color_hue(2)) + 
                    ggtitle(paste0(compn, '_', degt)) + 
                    theme(legend.position = 'bottom')

                print(fig, vp = viewport(x = 1/(nres*2) + (idx %% nres) * 1/nres, y = 1 - (1/(2*2) + (idx %/% nres) * 1/2), width = 1/nres, height = 1/2))
                idx <- idx + 1
            }
        }
        bouh <- dev.off()
    }


    ## build volcano plot for every comparison
    ctrl.volcdir <- paste0(figdir, "/Ctrl_AsyncVsSphase_comparisonVolcanoPlot")
    dir.create(ctrl.volcdir, showWarnings = F, recursive = F)
    
    volcano_plots(ctrl.resl, ctrl.volcdir)


    ## get summary about the whole set of DEGs
    if (TRUE) {
        ctrl.deg.summary <- deg.summary.fc(resl=ctrl.resl, expr.mtx=ctrl.vsd, fdr.thresh = fdr.thresh)

        ## save the summary
        ctrl.sfn <- paste0(resdir, '/ctrl_AsyncVsSphase_DEG_VSTexpr_dendrogram_summary.RDS')
        saveRDS(object = ctrl.deg.summary, file = ctrl.sfn)
    }
}



#####################################################
#### analysis of S-phase situation compared to asynchronous (~average) situation under effect of siASF1
#####################################################
if (TRUE) {
    cat.total.log <- desdf$RNA == "total"

    cat.spl_label_df <- data.frame(
        "sample" = c("D356T1", "D356T2", "D356T3", "D356T4", "D356T5", "D356T6", "D356T7", "D356T8", "D356T9", "D356T10", "D356T11", "D356T12", "D356T13", "D356T14", "D356T15", "D356T16"),
        "spl_label" = c("ctrl.0h.1", "ctrl.0h.2", "ctrl.2h.1", "ctrl.2h.2", "ctrl.5h.1", "ctrl.5h.2", "siASF1.0h.1", "siASF1.0h.2", "siASF1.2h.1", "siASF1.2h.2", "siASF1.5h.1", "siASF1.5h.2", "ctrl.async.1", "ctrl.async.2", "siASF1.async.1", "siASF1.async.2")
    )


    ## DESeq2 setting
    cat.dds <- DESeqDataSetFromMatrix(countData = mcnt[,which(cat.total.log) + 1],
                                colData = tibble(desdf[cat.total.log,]) %>% mutate(Sphase = ifelse(is.na(Sphase), 'none', as.character(Sphase))),
                                design = ~ condition + Sphase + condition:Sphase)
    cat.dds@colData@listData$Sphase <- factor(cat.dds@colData@listData$Sphase, levels=c('none', '0', '2', '5'))

    mcols(cat.dds) <- DataFrame(mcols(cat.dds), data.frame(Geneid = mcnt["Geneid"]))

    ## independent filtering
    keep <- rowSums(counts(cat.dds)) >= 10
    cat.dds <- cat.dds[keep,]


    ## DESeq processing
        # ???
    # cat.dds <- estimateSizeFactors(cat.dds)
    # cat.dds@colData@listData$sizeFactor <- cat.dds@colData@listData$sizeFactor * spi.fract.min.rel[cat.dds@colData@listData$sample]
    # cat.dds@colData@listData$sizeFactor <- spi.seqdep.nf[cat.dds@colData@listData$sample]
    # cat.dds <- estimateDispersions(cat.dds)
    # cat.dds <- nbinomWaldTest(cat.dds)

    cat.dds <- DESeq(cat.dds)
    cat.desq.outf <- paste0(resdir_deseq, '/', 'CtrlVsSiASF1_while_AsyncVsSphase_effect_DESeq2Results.RDS')
    saveRDS(cat.dds, cat.desq.outf)

    ## transforming data (VST)
    cat.vsd <- vst(cat.dds, blind=FALSE)


    ## distance heatmap/hclust of the samples
    sampleDists <- dist(t(assay(cat.vsd))) # default: method = "euclidean"
    sampleDistMatrix <- as.matrix(sampleDists)

    rownames(sampleDistMatrix) <- paste(cat.vsd$condition, cat.vsd$Sphase, sep=':')
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    figname <- paste0(figdir, "/CtrlVsSiASF1_while_AsyncVsSphase_effect_sample_correlation_hclust.pdf"); message(figname)
    pdf(figname)
    print(
        pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors)
    )
    bouh <- dev.off()


    ## exploring dataset structure
    figname <- paste0(figdir, "/CtrlVsSiASF1_AsyncVsSphase_effect_sample_PCA.pdf"); message(figname)
    pdf(figname)
    print(
        plotPCA(cat.vsd, intgroup=c("condition", "Sphase"))
    )
    bouh <- dev.off()


    ##  ## Explore histone gene expression along replication time points +/- ASF1
    if (TRUE) {
        hm_gene_vst_expression(vst.obj = cat.vsd, gene.sel.gff = histgr, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/histone_gene_expression"), spl_ref = "D356T1")
    }

    ##  ## for Shweta's article
    if (TRUE) {
        hist.cat.vsd <- cat.vsd[mcols(cat.vsd)$Geneid %in% histgr$gene_id,]

        sh.cat.vsd <- hist.cat.vsd[, paste0('D356T', c(15:16, 3:12))]

        ## treat for samples of synchronized cells
        #   # replicate 1
        assay(sh.cat.vsd[, paste0('D356T', c(3, 5, 7, 9, 11))]) <- apply(
            assay(sh.cat.vsd[, paste0('D356T', c(3, 5, 7, 9, 11))]),
            2,
            function(xxx) {
                xxx - unlist(assay(hist.cat.vsd[, 'D356T1']))
            }
        )

        #   # replicate 2
        assay(sh.cat.vsd[, paste0('D356T', c(4, 6, 8, 10, 12))]) <- apply(
            assay(sh.cat.vsd[, paste0('D356T', c(4, 6, 8, 10, 12))]),
            2,
            function(xxx) {
                xxx - unlist(assay(hist.cat.vsd[, 'D356T2']))
            }
        )

        # ## treat for samples of non-synchronized cells
        #   # replicate 1
        assay(sh.cat.vsd[, paste0('D356T', c(15))]) <- apply(
            assay(sh.cat.vsd[, paste0('D356T', c(15))]),
            2,
            function(xxx) {
                xxx - unlist(assay(hist.cat.vsd[, 'D356T13']))
            }
        )

        #   # replicate 2
        assay(sh.cat.vsd[, paste0('D356T', c(16))]) <- apply(
            assay(sh.cat.vsd[, paste0('D356T', c(16))]),
            2,
            function(xxx) {
                xxx - unlist(assay(hist.cat.vsd[, 'D356T14']))
            }
        )

            # ???
        hm_gene_vst_expression(vst.obj = sh.cat.vsd[,3:ncol(sh.cat.vsd)], gene.sel.gff = histgr_old, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/histone_gene_expression_synchro_ShwetaFig"), relative.values = TRUE)
        hm_gene_vst_expression(vst.obj = sh.cat.vsd[,1:2], gene.sel.gff = histgr_old, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/histone_gene_expression_async_ShwetaFig"), relative.values = TRUE, width.factor = 0.45)
    }

    if (TRUE) {
        fos.cat.vsd <- cat.vsd[mcols(cat.vsd)$Geneid %in% gene.gff$gene_id[(! is.na(gene.gff$gene_name)) & gene.gff$gene_name %in% c('FOS', 'GAPDH', 'H4C5', 'H1-4', 'H2BC7', 'H4C3')],]

        sf.cat.vsd <- fos.cat.vsd[, paste0('D356T', c(15:16, 3:12))]

        ## treat for samples of synchronized cells
        #   # replicate 1
        assay(sf.cat.vsd[, paste0('D356T', c(3, 5, 7, 9, 11))]) <- apply(
            assay(sf.cat.vsd[, paste0('D356T', c(3, 5, 7, 9, 11))]),
            2,
            function(xxx) {
                xxx - unlist(assay(fos.cat.vsd[, 'D356T1']))
            }
        )

        #   # replicate 2
        assay(sf.cat.vsd[, paste0('D356T', c(4, 6, 8, 10, 12))]) <- apply(
            assay(sf.cat.vsd[, paste0('D356T', c(4, 6, 8, 10, 12))]),
            2,
            function(xxx) {
                xxx - unlist(assay(fos.cat.vsd[, 'D356T2']))
            }
        )

        # ## treat for samples of non-synchronized cells
        #   # replicate 1
        assay(sf.cat.vsd[, paste0('D356T', c(15))]) <- apply(
            assay(sf.cat.vsd[, paste0('D356T', c(15))]),
            2,
            function(xxx) {
                xxx - unlist(assay(fos.cat.vsd[, 'D356T13']))
            }
        )

        #   # replicate 2
        assay(sf.cat.vsd[, paste0('D356T', c(16))]) <- apply(
            assay(sf.cat.vsd[, paste0('D356T', c(16))]),
            2,
            function(xxx) {
                xxx - unlist(assay(fos.cat.vsd[, 'D356T14']))
            }
        )

        disgr <- c(histgr_old[histgr_old$gene_name %in% c('HIST1H2BF', 'HIST1H4C', 'HIST1H1E', 'HIST1H4E')], gene.gff[(! is.na(gene.gff$gene_name)) & gene.gff$gene_name %in% c('FOS', 'GAPDH')])
        hm_gene_vst_expression(vst.obj = sf.cat.vsd[,3:ncol(sf.cat.vsd)], gene.sel.gff = disgr, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/smFISH_gene_expression_synchro_ShwetaFig"), relative.values = TRUE, height.factor = 0.13)
        hm_gene_vst_expression(vst.obj = sf.cat.vsd[,1:2], gene.sel.gff = disgr, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/smFISH_gene_expression_async_ShwetaFig"), relative.values = TRUE, width.factor = 0.45, height.factor = 0.13)
    }


    ##  ## Explore replicative histone generation paths
    if (TRUE) {
        hm_gene_vst_expression(vst.obj = cat.vsd, gene.sel.gff = replHistPathgr, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/replicative_histone_generation_gene_expression"), spl_ref = "D356T1")
        # hm_gene_vst_expression(vst.obj = counts(cat.dds), gene.sel.gff = replHistPathgr, spl_label_df = cat.spl_label_df, fig_base = paste0(figdir, "/replicative_histone_generation_gene_expression"), spl_ref = "D356T1", mat.gene.id = mcols(cat.dds)$Geneid)
    }


    # resultsNames(cat.dds)
    ## among the genes changing in expression between asynchronous and Sphase time point according to ASF1 activity, which genes show a variation in this behaviors across the time points
    cat.resl <- list(
        "AV0.CVsiA_vs_AV2.CVsiA" = results(cat.dds, contrast=list("conditionsiAsf1.Sphase0", "conditionsiAsf1.Sphase2")),
        "AV0.CVsiA_vs_AV5.CVsiA" = results(cat.dds, contrast=list("conditionsiAsf1.Sphase0", "conditionsiAsf1.Sphase5")),
        "AV2.CVsiA_vs_AV5.CVsiA" = results(cat.dds, contrast=list("conditionsiAsf1.Sphase2", "conditionsiAsf1.Sphase5"))
    )

    ## build volcano plot for every comparison
    cat.volcdir <- paste0(figdir, "/CtrlVsSiASF1_AsyncVsSphase_comparisonVolcanoPlot")
    dir.create(cat.volcdir, showWarnings = F, recursive = F)
    
    volcano_plots(cat.resl, cat.volcdir)


    ##  # look at the expression of the DEGs by comparison
    fdr.thresh <- 0.01
    l2fc.thresh <- 1.5 # only used for this chunck
    if (TRUE) {
        figname <- paste0(figdir, '/', 'DEG_expression_CtrlVsSiASF1_AsyncVsSphase_by_comparison_boxplot.pdf')
        pdf(figname, width = 1 * 7, height = 7)
        nres = length(cat.resl)
        grid.newpage()
        idx <- 0
        for (degt in c("up", "down")) {
            for (compn in names(cat.resl)) {
                dtab <- cat.vsd[(!is.na(cat.resl[[compn]]$padj)) & cat.resl[[compn]]$padj <= fdr.thresh & abs(cat.resl[[compn]]$log2FoldChange) >= l2fc.thresh & list("up" = cat.resl[[compn]]$log2FoldChange > 0, "down" = cat.resl[[compn]]$log2FoldChange < 0)[[degt]]]
                gdata <- tibble(gather(as.data.frame(assay(dtab)), sample, expression))
                gdata <- inner_join(gdata, desdf, by = "sample")
                gdata$sample_ <- paste0(gdata$condition, '_', gdata$Sphase)
                gdata$spl_name <- paste0(gdata$sample_, '_', gdata$replicate)
                gdata$spl_name <- factor(gdata$spl_name, levels = unique(gdata$spl_name), order = T)
                
                fig <- ggplot(gdata, aes(x = Sphase, y = expression, fill = condition, group = spl_name)) + 
                    geom_boxplot() + 
                    scale_fill_manual(values = gg_color_hue(2)) + 
                    ggtitle(paste0(compn, '_', degt)) + 
                    theme(legend.position = 'bottom')

                print(fig, vp = viewport(x = 1/(nres*2) + (idx %% nres) * 1/nres, y = 1 - (1/(2*2) + (idx %/% nres) * 1/2), width = 1/nres, height = 1/2))
                idx <- idx + 1
            }
        }
        bouh <- dev.off()
    }


    ## get summary about the whole set of DEGs
    if (TRUE) {
        cat.deg.summary <- deg.summary.fc(resl=cat.resl, expr.mtx=cat.vsd, fdr.thresh = fdr.thresh)

        ## save the summary
        cat.sfn <- paste0(resdir, '/ctrlVsSiASF1_AsyncVsSphase_DEG_VSTexpr_dendrogram_summary.RDS')
        saveRDS(object = cat.deg.summary, file = cat.sfn)

        ## recover results in variables
        cat.degil = cat.deg.summary[["res.deg.only"]]
        cat.mdegil = cat.deg.summary[["deg.log"]]
        cat.degvsd = cat.deg.summary[["deg.expr.mtx"]]
        cat.degvsd.center = cat.deg.summary[["center.expr.mtx"]]
        cat.deg.hclust = cat.deg.summary[["hclust"]]
        cat.deg.dendro = cat.deg.summary[["dendrogram"]]
    }


    ##  # plot gene expression of DEGs in every sample, with the Hclust clusters indicated
    if (TRUE) {
        nbr_clust_vec <- c(2, 4, 8, 10, 12, 16, 20, 24, 28, 32)
        nbr_clust_vec <- c(28, 29, 30, 31, 32)

        gdata <- tibble(gather(cbind(as.data.frame(cat.degvsd.center), "Geneid" = mcols(cat.dds)$Geneid[cat.mdegil]), sample, centered.expr, -Geneid))
        gdata$Geneid <- factor(x = gdata$Geneid,
            levels = mcols(cat.dds)$Geneid[cat.mdegil][order.dendrogram(cat.deg.dendro)], 
            ordered = TRUE)

        gdata$sample <- factor(x = gdata$sample,
            levels = unique(gdata$sample), 
            ordered = TRUE)

        deg.heatmap.plot <- ggplot(data = gdata, 
            aes(x = Geneid, y = sample)) +
            geom_tile(aes(fill = centered.expr)) +
            # scale_fill_viridis_c() +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
            theme(axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                # axis.text.y = element_blank(),
                # axis.title.y = element_blank(),
                # axis.ticks.y = element_blank(),
                legend.position = "left")



        cat.hcll <- list()
        for (kkk in nbr_clust_vec) {
            cat.deg.grp <- cutree(tree = cat.deg.hclust, k = kkk)
            dgdata <- tibble(data.frame("Geneid" = mcols(cat.vsd)$Geneid[cat.mdegil], "group" = factor(cat.deg.grp), "dendro.order" = factor(order(order.dendrogram(cat.deg.dendro)))))
            cat.hcll[[as.character(kkk)]] <- ggplot(dgdata, aes(x = dendro.order, y = 1, fill = group)) + 
                                geom_tile() +
                                scale_fill_manual(values = sample(rainbow(kkk), size = kkk)) +
                                theme(
                                    axis.text.x = element_blank(),
                                    axis.title.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.title.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    legend.position = 'none'
                                )
        }

        if (TRUE) {
            figname <- paste0(figdir, '/', 'ctrlVsSiASF1_AsyncVsSphase_DEG_centeredExpr_hclustLab_heatmap.pdf'); message(figname)
            pdf(figname, height = 7 * 1, width = 7 * 4)
            grid.newpage()
            print(deg.heatmap.plot, vp = viewport(x = 0.5, y = 0.75, width = 1, height = 0.5))
            for (iii in seq_len(length(nbr_clust_vec))) {
                kkk <- nbr_clust_vec[iii]
                print(
                    cat.hcll[[as.character(kkk)]],
                    vp = viewport(x = 0.535, y = 1/(2^2 * length(nbr_clust_vec)) + ((iii - 1) * (1 /(2 * length(nbr_clust_vec)))), width = 0.930, height = 1/(2 * length(nbr_clust_vec)))
                )
            }
            bouh <- dev.off()
        }
    }


    #### GO Term enrichment on siASF1 DEG clusters in asynchronous cell culture
    #### by ORA (on each cluster of DEGs)
    nclust <- 30
    if (TRUE) {
        ## compute enrichment by ORA
        cat.gedir <- paste0(resdir, '/geneEnrichAnalysis_CtrlVsSiASF1_AsyncVsSphase')
        cat.oral <- ora_go_enrichment(cat.deg.hclust, degids = mcols(cat.degvsd)$Geneid, allids = mcols(cat.vsd)$Geneid, nclust = nclust, gedir = cat.gedir)

        ## recover results in variables
        cat.deg.grp = cat.oral[["deg.cluster.id"]]
        cat.erchl = cat.oral[["go.enrichment"]]

        ## save enrichment analysis (because it takes a lot of time to compute)
        cat.ofn <- paste0(resdir, '/enrichment_ORA_15HierarchicalClustersOfDEGs_CtrlVsSiASF1_AsyncVsSphase_variousPathways.RDS')
        saveRDS(object = cat.erchl, file = cat.ofn)

        ## path for ORA plots
        cat.enrdir <- paste0(figdir, "/GOenrichment_CtrlVsSiASF1_AsyncVsSphase")
        dir.create(cat.enrdir, showWarnings = F, recursive = F)

        sclist <- list(
            "geneontology_Biological_Process_noRedundant" = "BP",
            "geneontology_Cellular_Component_noRedundant" = "CC",
            "pathway_KEGG" = "KEGG",
            "pathway_Wikipathway" = "WP"
        )

        ## build the plots
        for (cii in names(cat.erchl)) {
            cat.erch <- cat.erchl[[cii]]

            for (gset in names(cat.erch)) {
                # gset <- names(cat.erch)[1]

                gerch <- tibble(cat.erch[[gset]])
                gerch$description <- factor(gerch$description, levels=unique(gerch$description[order(gerch$enrichmentRatio)]))
                sgerch <- gerch[gerch$FDR < 0.05,]

                if (nrow(sgerch) > 0) {
                    fig <- ggplot(
                        sgerch
                    ) +
                        geom_segment(aes(x=description, y=0, xend=description, yend=enrichmentRatio)) +
                        geom_point(aes(x=description, y=enrichmentRatio), size=5) +
                        geom_hline(yintercept=1.5) +
                        xlab("GO term") +
                        ylab("Enrichment Ratio") +
                        coord_flip() +
                        theme(text=element_text(size=20))

                    figname <- paste0(cat.enrdir, "/cluster", cii, "_", sclist[[gset]], ".pdf"); message(figname)

                    pdf(figname, width = 7 * 1.5)
                    print(fig)
                    bouh <- dev.off()
                }
            }
        }

        ## build tables with genes making the GO enrichment by cluster
        cat.resdf <- go_term_enrich_genes(erchl = cat.erchl, erch.fdr = 0.05, go.db = "geneontology_Biological_Process_noRedundant", mart = mart)

        cat.resdffn <- paste0(resdir, "/GOenrichment_CtrlVsSiASF1_AsyncVsSphase_geneTable.tsv")
        write.table(cat.resdf, cat.resdffn, quote=F, row.names=F, col.names=T, sep='\t')


    }
}



#####################################################
#### Looking for potential genes impaired in their expression homeostasis
#####################################################
if (TRUE) {
    #### Whole set of synchrosized cells
    ehm.dds.vsd <- dds.analysis(mcnt[,2:13], desdf[1:12,], Geneid = mcnt["Geneid"], design = ~ condition + Sphase + Sphase:condition)
    ehm.dds = ehm.dds.vsd[["dds"]]
    ehm.vsd = ehm.dds.vsd[["vsd"]]
    ehm.spl_label_df <- data.frame(
        "sample" = c("D356T1", "D356T2", "D356T3", "D356T4", "D356T5", "D356T6", "D356T7", "D356T8", "D356T9", "D356T10", "D356T11", "D356T12"),
        "spl_label" = c("ctrl.0h.1", "ctrl.0h.2", "ctrl.2h.1", "ctrl.2h.2", "ctrl.5h.1", "ctrl.5h.2", "siASF1.0h.1", "siASF1.0h.2", "siASF1.2h.1", "siASF1.2h.2", "siASF1.5h.1", "siASF1.5h.2")
    )

    # resultsNames(dds)
    ehm.resl <- list(
        "Asf1VsCtrl" = results(ehm.dds, contrast=c("condition", "siAsf1", "siCtrl")), # main effect of siAsf1 vs siCtrl on 0h of S-phase
        "A2VsA0" = results(ehm.dds, contrast=c("Sphase", "2", "0")), # main effect of S-phase at 2h vs 0h in condition of siCtrl
        "5Vs0" = results(ehm.dds, contrast=c("Sphase", "5", "0")), # main effect of S-phase at 5h vs 0h in condition of siCtrl
        "5Vs2" = results(ehm.dds, contrast=c("Sphase", "5", "2")), # main effect of S-phase at 5h vs 2h in condition of siCtrl
        "5From2Vs0" = results(ehm.dds, contrast=list(c("Sphase_5_vs_0", "Sphase_2_vs_0"))), # additional effect of S-phase at 5h to effet of 2h vs 0h in condition of siCtrl
        "A2VsC0" = results(ehm.dds, name="conditionsiAsf1.Sphase2"), # test if the effec of S-phase at 2h vs 0h is different in condition siAsf1 compared to siCtrl
        "A5VsC0" = results(ehm.dds, name="conditionsiAsf1.Sphase5"), # test if the effec of S-phase at 5h vs 0h is different in condition siAsf1 compared to siCtrl
        "A5FromA2VsC0" = results(ehm.dds, contrast=list(c("conditionsiAsf1.Sphase5", "conditionsiAsf1.Sphase2"))) # variation between 5h compared to 2h vs 0h of S-phase in the additional effect by siAsf1 vs siCtrl.
    )


    #### siCtrl set of synchronized cells
    con.dds.vsd <- dds.analysis(mcnt[,2:7], desdf[1:6,], Geneid = mcnt["Geneid"], design = ~ Sphase)
    con.dds = con.dds.vsd[["dds"]]
    con.vsd = con.dds.vsd[["vsd"]]
    con.spl_label_df <- data.frame(
        "sample" = c("D356T7", "D356T8", "D356T9", "D356T10", "D356T11", "D356T12"),
        "spl_label" = c("ctrl.0h.1", "ctrl.0h.2", "ctrl.2h.1", "ctrl.2h.2", "ctrl.5h.1", "ctrl.5h.2")
    )

    # resultsNames(dds)
    con.resl <- list(
        "2Vs0" = results(con.dds, contrast=c("Sphase", "2", "0")),
        "5Vs2" = results(con.dds, contrast=c("Sphase", "5", "2")),
        "5Vs0" = results(con.dds, contrast=c("Sphase", "5", "0"))
    )


    #### siAsf1 set of synchronized cells
    aon.dds.vsd <- dds.analysis(mcnt[,8:13], desdf[7:12,], Geneid = mcnt["Geneid"], design = ~ Sphase)
    aon.dds = aon.dds.vsd[["dds"]]
    aon.vsd = aon.dds.vsd[["vsd"]]
    aon.spl_label_df <- data.frame(
        "sample" = c("D356T1", "D356T2", "D356T3", "D356T4", "D356T5", "D356T6"),
        "spl_label" = c("siASF1.0h.1", "siASF1.0h.2", "siASF1.2h.1", "siASF1.2h.2", "siASF1.5h.1", "siASF1.5h.2")
    )

    # resultsNames(dds)
    aon.resl <- list(
        "2Vs0" = results(aon.dds, contrast=c("Sphase", "5", "0")),
        "5Vs2" = results(aon.dds, contrast=c("Sphase", "5", "2")),
        "5Vs0" = results(aon.dds, contrast=c("Sphase", "5", "0"))
    )


    #### Look for genes with overexpression across replication time points in siAsf1 condition while staying at same expression in siCtrl condition
    if (TRUE) {
        aon.deg.up.log <- lapply(aon.resl, function(xxx) {
            xxx$log2FoldChange > 0 & xxx$padj < 0.01
        })
        aon.deg.up.log <- apply(
            do.call(cbind, aon.deg.up.log),
            1,
            all,
            na.rm = F
        )
        aon.deg.up.log[is.na(aon.deg.up.log)] <- FALSE
    }

    if (TRUE) {
        con.ceg.log <- lapply(con.resl, function(xxx) {
            xxx$log2FoldChange < 0.5 & xxx$padj > 0.15
        })
        con.ceg.log <- apply(
            do.call(cbind, con.ceg.log),
            1,
            all,
            na.rm = F
        )
        con.ceg.log[is.na(con.ceg.log)] <- FALSE
    }

    # selection of the genes of interest among the DEGs in siASF1 and in siCtrl samples (up-reg in siASF1 and stable in siContr)
    ehm.deg.log <- aon.deg.up.log[mcols(aon.dds)$Geneid %in% mcols(con.dds)$Geneid] & con.ceg.log[mcols(con.dds)$Geneid %in% mcols(aon.dds)$Geneid]
    # get select of the same genes in the whole set of genes
    ehm.log <- mcols(ehm.dds)$Geneid %in% mcols(aon.dds)$Geneid[
        mcols(aon.dds)$Geneid %in% mcols(con.dds)$Geneid
        ] & mcols(ehm.dds)$Geneid %in% mcols(con.dds)$Geneid[
                    mcols(con.dds)$Geneid %in% mcols(aon.dds)$Geneid
                    ]


    # expression in heatmap
    if (TRUE) {
        ehm.toutp <- paste0(figdir, "/siCtrlConstant_siASF1UpThroughReplication_gene_expression")
        dir.create(ehm.toutp, showWarnings = F, recursive = F)
        gene.sel.gr <- gene.gff[(!is.na(gene.gff$gene_id)) & gene.gff$gene_id %in% mcols(ehm.dds)$Geneid[ehm.log][ehm.deg.log]]
        hm_gene_vst_expression(vst.obj = ehm.vsd, gene.sel.gff = gene.sel.gr, spl_label_df = ehm.spl_label_df, fig_base = paste0(ehm.toutp, "/CtrlAndSiASF1_Sphase_gene_expression"), spl_ref = 'D356T1', height.factor = length(gene.sel.gr) / 40)
    }


    #### cluster genes of interest based on "optimal transport" (~ Wasserstein/Kantorovitch) distance
    if (TRUE) {
        wrk.ehm.vsd <- ehm.vsd[mcols(ehm.vsd)$Geneid %in% gene.sel.gr$gene_id]
        assay(wrk.ehm.vsd) <- assay(wrk.ehm.vsd) - rowMeans(assay(wrk.ehm.vsd))  # / rowSums(assay(wrk.ehm.vsd))))

        dstr.list <- data.frame(t(assay(wrk.ehm.vsd)))
        gene.name.vec <- left_join(data.frame(mcols(wrk.ehm.vsd)), data.frame(mcols(gene.sel.gr)) %>% mutate(Geneid = gene_id), by = 'Geneid')[['gene_name']]
        colnames(dstr.list) <- gene.name.vec
        dstr.list <- as.list(dstr.list)
        
        wass.dist.mat <- function(dstr.list) {
            dmtx <- matrix(0, length(dstr.list), length(dstr.list))
            
            for (xii in 1:(length(dstr.list) - 1)) {
                for (yii in (xii + 1):(length(dstr.list))) {
                    dmtx[yii, xii] <- wasserstein1d(dstr.list[[xii]], dstr.list[[yii]])
                }
            }

            return(as.dist(dmtx))
        }

        wdmtx <- wass.dist.mat(dstr.list)
        
        ehm.wh <- hclust(wdmtx)
        ehm.wh.dendro  <- as.dendrogram(ehm.wh)
        plot(ehm.wh.dendro)

        dev.new()
        gdata <- reshape2::melt(assay(wrk.ehm.vsd))
        colnames(gdata) <- c("gene", "sample", "vst.expr")
        gdata$gene <- gene.name.vec

        ggplot(
            gdata %>% mutate(sample = factor(sample, levels = rev(unique(sample))), gene = factor(gene, level = unique(gene)[order.dendrogram(ehm.wh.dendro)])),
            aes(
                x = gene,
                y = sample,
                fill = vst.expr
            )
        ) +
            geom_tile() +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red")
    }
}



#####################################################
#####################################################
#### Analyse by DESeq2 the sample of nascent RNAs
#####################################################
    #### Whole set of nascent RNA samples
    nas.nascent.log <- desdf$RNA == "nascent"

    nas.dds.vsd <- dds.analysis(mcnt[,which(nas.nascent.log) + 1], desdf[nas.nascent.log,], Geneid = mcnt["Geneid"], design = ~ condition)
    nas.dds = nas.dds.vsd[["dds"]]
    nas.vsd = nas.dds.vsd[["vsd"]]
    nas.spl_label_df <- data.frame(
        "sample" = c("D356T17", "D356T18", "D356T19", "D356T20"),
        "spl_label" = c("ctrl.nasc.1", "siASF1.nasc.1", "ctrl.nasc.2", "siASF1.nasc.2")
    )

    ## Saving DESeq2 processing
    nas.desq.outf <- paste0(resdir_deseq, '/', 'CtrlVsSiASF1_inNascentRNAs_DESeq2Results.RDS')
    saveRDS(nas.dds, nas.desq.outf)


    ##  ## Explore histone gene expression +/- ASF1
    if (TRUE) {
        hm_gene_vst_expression(vst.obj = nas.vsd, gene.sel.gff = histgr, spl_label_df = nas.spl_label_df, fig_base = paste0(figdir, "/histone_gene_expression_nascent"), spl_ref = NULL, width.factor = 0.5)
    }

    ##  ## for Shweta's article
    if (TRUE) {
        hist.nas.vsd <- nas.vsd[mcols(nas.vsd)$Geneid %in% histgr$gene_id,]

        sh.nas.vsd <- hist.nas.vsd[, paste0('D356T', c(18, 20))]

        ## normalize siAsf1 samples by siCtrl samples
        #   # replicate 1
        assay(sh.nas.vsd[, 'D356T18']) <- apply(
            assay(sh.nas.vsd[, 'D356T18']),
            2,
            function(xxx) {
                xxx - unlist(assay(hist.nas.vsd[, 'D356T17']))
            }
        )

        #   # replicate 2
        assay(sh.nas.vsd[, 'D356T20']) <- apply(
            assay(sh.nas.vsd[, 'D356T20']),
            2,
            function(xxx) {
                xxx - unlist(assay(hist.nas.vsd[, 'D356T19']))
            }
        )

        hm_gene_vst_expression(vst.obj = sh.nas.vsd, gene.sel.gff = histgr_old, spl_label_df = nas.spl_label_df, fig_base = paste0(figdir, "/histone_gene_expression_nascent_ShwetaFig"), relative.values = TRUE, width.factor = 0.45)
    }

    if (TRUE) {
        fos.nas.vsd <- nas.vsd[mcols(nas.vsd)$Geneid %in% gene.gff$gene_id[(! is.na(gene.gff$gene_name)) & gene.gff$gene_name %in% c('FOS', 'GAPDH', 'H4C5', 'H1-4', 'H2BC7', 'H4C3')],]

        sf.nas.vsd <- fos.nas.vsd[, paste0('D356T', c(18, 20))]

        ## normalize siAsf1 samples by siCtrl samples
        #   # replicate 1
        assay(sf.nas.vsd[, 'D356T18']) <- apply(
            assay(sf.nas.vsd[, 'D356T18']),
            2,
            function(xxx) {
                xxx - unlist(assay(fos.nas.vsd[, 'D356T17']))
            }
        )

        #   # replicate 2
        assay(sf.nas.vsd[, 'D356T20']) <- apply(
            assay(sf.nas.vsd[, 'D356T20']),
            2,
            function(xxx) {
                xxx - unlist(assay(fos.nas.vsd[, 'D356T19']))
            }
        )

        disgr <- c(histgr_old[histgr_old$gene_name %in% c('HIST1H2BF', 'HIST1H4C', 'HIST1H1E', 'HIST1H4E')], gene.gff[(! is.na(gene.gff$gene_name)) & gene.gff$gene_name %in% c('FOS', 'GAPDH')])
        hm_gene_vst_expression(vst.obj = sf.nas.vsd, gene.sel.gff = disgr, spl_label_df = nas.spl_label_df, fig_base = paste0(figdir, "/smFISH_gene_expression_nascent_ShwetaFig"), relative.values = TRUE, width.factor = 0.45, height.factor = 0.13)
    }

#####################################################

########

########


#####################################################
#####################################################
####
## Differential gene expression analysis of Schweta's project (+/- ASF1)

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

oldtheme <- theme_set(theme_bw())

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

if (TRUE) {
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
#### analysis of S-phase situation compared to asynchronous (~average) situation under effect of siASF1
#####################################################
if (TRUE) {
    cat.total.log <- desdf$RNA == "total"

    cat.spl_label_df <- data.frame(
        "sample" = c("D356T1", "D356T2", "D356T3", "D356T4", "D356T5", "D356T6", "D356T7", "D356T8", "D356T9", "D356T10", "D356T11", "D356T12", "D356T13", "D356T14", "D356T15", "D356T16"),
        "spl_label" = c("ctrl.0h.1", "ctrl.0h.2", "ctrl.2h.1", "ctrl.2h.2", "ctrl.5h.1", "ctrl.5h.2", "siASF1.0h.1", "siASF1.0h.2", "siASF1.2h.1", "siASF1.2h.2", "siASF1.5h.1", "siASF1.5h.2", "ctrl.async.1", "ctrl.async.2", "siASF1.async.1", "siASF1.async.2")
    )


    ## DESeq2 setting
    cat.dds <- DESeqDataSetFromMatrix(countData = mcnt[,which(cat.total.log) + 1],
                                colData = tibble(desdf[cat.total.log,]) %>% mutate(Sphase = ifelse(is.na(Sphase), 'none', as.character(Sphase))),
                                design = ~ condition + Sphase + condition:Sphase)
    cat.dds@colData@listData$Sphase <- factor(cat.dds@colData@listData$Sphase, levels=c('none', '0', '2', '5'))

    mcols(cat.dds) <- DataFrame(mcols(cat.dds), data.frame(Geneid = mcnt["Geneid"]))

    ## independent filtering
    keep <- rowSums(counts(cat.dds)) >= 10
    cat.dds <- cat.dds[keep,]


    ## DESeq processing
    cat.dds <- DESeq(cat.dds)
    cat.desq.outf <- paste0(resdir_deseq, '/', 'CtrlVsSiASF1_while_AsyncVsSphase_effect_DESeq2Results.RDS')
    saveRDS(cat.dds, cat.desq.outf)
}



#####################################################
#####################################################
#### Analyse by DESeq2 the sample of nascent RNAs
#####################################################
    #### Whole set of nascent RNA samples
    nas.nascent.log <- desdf$RNA == "nascent"

    nas.dds.vsd <- dds.analysis(mcnt[,which(nas.nascent.log) + 1], desdf[nas.nascent.log,], Geneid = mcnt["Geneid"], design = ~ condition)
    nas.dds = nas.dds.vsd[["dds"]]
    nas.vsd = nas.dds.vsd[["vsd"]]
    nas.spl_label_df <- data.frame(
        "sample" = c("D356T17", "D356T18", "D356T19", "D356T20"),
        "spl_label" = c("ctrl.nasc.1", "siASF1.nasc.1", "ctrl.nasc.2", "siASF1.nasc.2")
    )

    ## Saving DESeq2 processing
    nas.desq.outf <- paste0(resdir_deseq, '/', 'CtrlVsSiASF1_inNascentRNAs_DESeq2Results.RDS')
    saveRDS(nas.dds, nas.desq.outf)

#####################################################
