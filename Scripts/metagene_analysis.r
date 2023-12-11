
{
    test <- TRUE
}

if (TRUE) {
    require(rtracklayer)
    require(tidyverse)
    require(data.table)

    oldtheme <- theme_set(theme_bw())
    theme_txt <- theme(text = element_text(size = 20))
    theme_txt_xangle <- theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust=1))
}


if (TRUE) {
    source('./Scripts/general/_functions_in_R.r')

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


    bin_signal <- function(dwm, extra_res) {
        ## get the number of elements that fit in the resolution (every bin must have all the atomic values of the region it summarizes)
        nebin <- floor(length(dwm) / extra_res)

        ## reduce the vector to values allowing to completely fill every bin
        dwm <- dwm[1:(extra_res * nebin)]

        ## set which value to which bin
        colix <- floor((1:length(dwm) - 1) / extra_res) + 1
        colixn <- tapply(1:length(dwm), colix, identity)

        ## summarize values in every bin
        for (colixi in 1:nebin) {
            dwm[colixi] <- mean(dwm[colixn[[colixi]]])
        }
        dwm <- dwm[1:nebin]

        ## return results
        return(dwm)
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

    ## color palette for conditions
    treat_colpal <- c('siCtrl' = '#3174a1', 'siAsf1' = '#c03d3e')

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

    ## load curated list of annotation versions of replicative histone genes
    rplhgr <- rtracklayer::import.gff("./Data/ASF1_chaperone/annotations/replicative_histone_genes.gff")

    ## load annotations of Stem Loop (SL)
    slgr <- rtracklayer::import.gff("./Results/ASF1_chaperone/transcriptome_signal_in_SL_and_HDE/histone_stem_loop.gff")

    ## load list of genes involved in replicative histone biogenesis
    replHistPathgr <- rtracklayer::import.gff('./Data/path_replicative_histone_synthesis.gtf')

    ## load clustering of the genes
    albcl.p <- './Data/ASF1_chaperone/Alberto_Gatto_results/results/clusters_6/gene_clusters_sync.tsv'
    albcl.t <- read.table(albcl.p, h=T, sep='\t')
    colnames(albcl.t)[1] <- 'gene_id'
    # head(albcl.t)
    gr2g <- albcl.t$name[albcl.t$cluster == 2]
    gr2g <- hist_name_oldnew$new[hist_name_oldnew$old %in% gr2g]


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
        "replicate" = factor(as.character(c(rep(1:2, 8), rep(1:2, each=2))))
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

    ## load spike and sequencing depth normalization factors
    spi.seqdep.nf.path <- file.path(resdir, "spike_and_seqDepth_norm_factor.tsv")
    spi.seqdep.nf <- read.table(spi.seqdep.nf.path, h=T, sep='\t')

    ## load spike normalization factors
    spi.nf.path <- file.path(resdir, "spike_norm_factor.tsv")
    spi.nf <- read.table(spi.nf.path, h=T, sep='\t')
}


#####################################################
#### metagene
#####################################################
if (TRUE) {

    #### set the list of annotations for metagene
    {
        ## metagene parameters
        {
            upstream_length <- 50
            downstream_length <- 500
            nbin <- 20

            mtg_width <- 500
            extra_res <- 10
        }

        ## restrict to Group 2 genes
        {
            srplhgr <- rplhgr[rplhgr$gene_name %in% gr2g]
            sslgr <- slgr[slgr$gene_id %in% gr2g]
        }

        ## remove neighbour convergent genes H2AC21 and H2AC20
        {
            srplhgr <- srplhgr[!(srplhgr$gene_name %in% c('H2AC21', 'H2AC20'))]
            sslgr <- sslgr[!(sslgr$gene_id %in% c('H2AC21', 'H2AC20'))]
        }

        ## select relplicative histone genes with stem loop
        {
            srplhgr <- srplhgr[srplhgr$gene_name %in% sslgr$gene_id]
        }


        ## sort replicative histone genes by the stem loops
        {
            srplhgr$sl_ordid <- left_join(
                data.frame(gene_name=srplhgr$gene_name),
                data.frame(gene_name=sslgr$gene_id, ordid=seq(1, length(sslgr))),
                by = 'gene_name'
            )$ordid

            # srplhgr <- rplhgr[!is.na(rplhgr$sl_ordid)]
            srplhgr <- srplhgr[order(srplhgr$sl_ordid)]

        }


        ## initiate GRanges of metagene annotations from the SL coordinates
        {
            anngr <- GRanges(
                seqnames = seqnames(sslgr),
                ranges = IRanges(start=start(sslgr), end=end(sslgr)),
                strand=strand(sslgr),
                seqinfo=seqinfo(sslgr),
                seqlengths=seqlengths(sslgr),
                gene_name = sslgr$gene_id
            )
        }

        ## update the annotation borders according the strand
        {
            ## get logicals by strand
            {
                log_strand_pos <- as.vector((strand(anngr) == '+'))
                log_strand_neg <- !log_strand_pos
            }

            ## extend to TSS of the gene
            {
                start(anngr)[log_strand_pos] <- start(srplhgr)[log_strand_pos]
                end(anngr)[log_strand_neg] <- end(srplhgr)[log_strand_neg]
                # end(anngr) - start(anngr)
            }

            ## extend upstream anddownstream
            {
                start(anngr)[log_strand_pos] <- start(anngr)[log_strand_pos] - upstream_length
                end(anngr)[log_strand_neg] <- end(anngr)[log_strand_neg] + upstream_length

                end(anngr)[log_strand_pos] <- end(anngr)[log_strand_pos] + downstream_length
                start(anngr)[log_strand_neg] <- start(anngr)[log_strand_neg] - downstream_length
            }
        }
    }


    #### get the signal per annotation
    if (FALSE)
    {
        ## compute the metaplot signals
        cust_metagene_signal_fc <- function(anngr, desdf, upstream_length, downstream_length, extra_res, nbin, mtg_width, spi.seqdep.nf)
        {
            tots <- list()
            covdfl <- list()
            for (spln in desdf$sample)
            {
                # spln <- "D356T17"

                ## load the coverage track (no 0 values)
                bwp <- file.path("Data", "ASF1_chaperone", "RNA-Seq", "coverage", paste0(spln, ".sorted.markdup.bw"))
                message(bwp)
                # file.exists(bwp)
                covgr <- rtracklayer::import.bw(bwp)

                ## recover signal along the metagene annotations
                {
                    ## find coverage blocks per replicative histone gene
                    ovl <- findOverlaps(covgr, anngr)
                    ovlt <- tapply(queryHits(ovl), subjectHits(ovl), identity)

                    ## compute the coverage track from the blocks
                    covl <- list()
                    tots[[spln]] <- list()
                    for (anni in 1:length(anngr)) {
                        # anni <- 1

                        ## extract the coverage blocks
                        scovgr <- covgr[ovlt[[anni]]]
                        start(scovgr)[1] <- start(anngr[anni])
                        end(scovgr)[length(scovgr)] <- end(anngr[anni])
                        scovgr$ordi <- as.numeric(1:length(scovgr))

                        ## identify the missing 0-value blocks
                        {
                            endv <- end(scovgr)
                            startv <- start(scovgr)
                            spacev <- tail(startv, n=-1) - head(endv, n=-1) - 1
                            log_zbl <- spacev != 0
                            if (any(log_zbl)) {
                                id_zbl <- which(log_zbl)

                                zgr <- scovgr[id_zbl]
                                zgr$score <- 0
                                zgr$ordi <- zgr$ordi + 0.5
                                end(zgr) <- start(scovgr)[id_zbl + 1] - 1
                                start(zgr) <- end(scovgr)[id_zbl] + 1

                                scovgr <- c(scovgr, zgr)
                                scovgr <- scovgr[order(scovgr$ordi)]
                            }
                        }


                        ## build coverage vector 
                        covv <- rep(scovgr$score, width(scovgr))
                        # print(c(tail(end(scovgr), n=1) - start(scovgr)[1] + 1, length(covv), width(anngr[anni])))
                        # if (tail(end(scovgr), n=1) - start(scovgr)[1] + 1 != length(covv)) {
                        #     message('longueurs différentes')
                        #     print(anni)
                        # }

                        ## normalize by spike-in controls and sequencing depth
                        covv <- spi.seqdep.nf$norm_factor[spi.seqdep.nf$sample == spln] * covv

                        ## get total signal in metagene normalized by gene length (before transformation in metagene)
                        tots[[spln]][[anngr$gene_name[anni]]] <- mean(covv[(upstream_length + 1):(upstream_length + width(anngr[anni]) - upstream_length - downstream_length)])

                        ## adjust by transcription sens
                        if (as.vector(strand(anngr[anni])) == '-') {
                            covv <- rev(covv)
                        }


                        ## transform to metagene scale
                        {
                            ## partition the signal
                            {
                                upv <- head(covv, upstream_length)
                                mtgv <- covv[(upstream_length + 1):(length(covv) - downstream_length)]
                                dwv <- tail(covv, downstream_length)
                                # all(c(upv, mtgv, dwv) == covv)
                            }

                            ## summarize the metagene part
                            {
                                binsz <- length(mtgv) / nbin
                                brds <- floor(binsz * 1:nbin)
                                ubrds <- c(1, head(brds, n=-1) + 1)

                                mtgbv <- apply(data.frame(u=ubrds, d=brds), 1, function(xxx, vec) {
                                    out <- vec[xxx[1]:xxx[2]]
                                    out <- mean(out)

                                    return(out)
                                }, vec=mtgv)
                            }
                        }

                        ## record the metagene track
                        covl[[anngr$gene_name[anni]]] <- c(upv, mtgbv, dwv)

                    }

                    ## reshape in data.frame
                    {
                        covdf <- as.data.frame(do.call(rbind, covl))
                        remove(covl)

                        colnames(covdf) <- c(paste0('u', 1:upstream_length), paste0('m', 1:nbin), paste0('d', 1:downstream_length))
                        covdfl[[spln]] <- covdf
                        remove(covdf)
                    }
                }
            }


            ## normalize by total signal in the replicative histone gene in siCtrl condition
            {
                totsdf <- reshape2::melt(tots)
                colnames(totsdf) <- c('signal', 'gene_name', 'sample')
                # str(totsdf)
                stotsdf = totsdf
                totsdf$condition <- left_join(
                    totsdf['sample'],
                    desdf,
                    by = 'sample'
                )$condition
                totsdf <- pivot_wider(totsdf, names_from='gene_name', values_from='signal')

                ctots <- colSums(totsdf[totsdf$condition == "siCtrl", anngr$gene_name])
                ctots.log <- log10(ctots)
                ctots.log.avg <- mean(ctots.log)
                ctots.log.difftoavg <- ctots.log - ctots.log.avg
                sig.cds.nf <- 10^(-ctots.log.difftoavg)

                ncovdfl <- list()
                for (spln in names(covdfl)) {
                    covdf = covdfl[[spln]]

                    ncovdf <- as.data.frame(lapply(covdf, function(xxx, sig.cds.nf) {
                        out = xxx * sig.cds.nf
                        return(out)
                    }, sig.cds.nf))

                    ncovdfl[[spln]] <- ncovdf
                }
            }


            ## transform raw signal in metagene track
            mtgcovl <- list()
            for (spln in names(ncovdfl))
            {
                ncovdf = ncovdfl[[spln]]

                ## plot metagene signal per replicative histone
                if (FALSE) {
                    gdata <- pivot_longer(ncovdf, cols = colnames(ncovdf), values_to = 'signal', names_to = 'name')
                    gdata$gene <- rep(rownames(ncovdf), each=ncol(ncovdf))

                    ## coordinate
                    coordv <- seq(-upstream_length, -1, 1)
                    coordv <- c(coordv, tail(coordv, 1) + (1:nbin) * (500 / nbin))
                    coordv <- c(coordv, tail(coordv, 1) + (1:downstream_length))
                    gdata$coordv <- rep(coordv, nrow(ncovdf))

                    ## plot
                    fig <- ggplot(
                        gdata,
                        aes(
                            x = coordv,
                            y = signal,
                            col = gene
                        )
                    ) +
                        geom_line() + 
                        theme_txt
                    
                    figname <- file.path(figdir, paste0('test_', spln, '.pdf')); message(figname)
                    pdf(figname, width = 2 * 7)
                    print(fig)
                    bouh <- dev.off()
                }

                ## compute metacoverage
                mtgcov <- unlist(lapply(ncovdf, mean))

                ## summarize extra regions in bins
                {
                    ## downstream region
                    {
                        downstream_length_fin <- floor(downstream_length / extra_res)
                        upstream_length_fin <- floor(upstream_length / extra_res)

                        dwm <- bin_signal(mtgcov[grepl('^d', names(ncovdf))], extra_res)
                        upm <- rev(bin_signal(rev(mtgcov[grepl('^u', names(ncovdf))]), extra_res))

                        ## adjust the names of the values
                        names(dwm) <- paste0('d', ((0:(downstream_length_fin - 1)) + 1) * extra_res)
                        names(upm) <- paste0('u', - ((0:(upstream_length_fin - 1)) + 1) * extra_res)

                        ## combine the regions
                        mtgcov <- c(upm, mtgcov[grepl('^m', names(mtgcov))], dwm)
                    }
                }

                ## reshape in dataframe
                mtgcov <- named_vec_to_df(mtgcov)
                remove(ncovdf)

                ## set metagene coordinates
                {
                    coordv <- (-upstream_length_fin:-1) * extra_res
                    coordv <- c(coordv, 0:(nbin - 1) * (mtg_width / nbin))
                    coordv <- c(coordv, mtg_width + ((1:downstream_length_fin) - 1) * extra_res)
                    mtgcov$coord <- coordv
                    remove(coordv)
                }

                ## record metagene signal for the sample
                mtgcovl[[spln]] <- mtgcov
                remove(mtgcov)
            }


            ## normalize by replicate effect
            {
                ## get total signal by sample in metagene
                rpsigl <- lapply(mtgcovl, function(xxx) {
                    totsig <- sum(xxx$value[grepl('^m', xxx$name)])

                    return(totsig)
                })

                rpsigdf <- reshape2::melt(rpsigl)
                colnames(rpsigdf) <- c('value', 'sample')

                rpsigdf$replicate <- left_join(
                    rpsigdf['sample'],
                    desdf,
                    by = 'sample'
                )$replicate


                ## get normalization factor per replicate
                avg_repsig <- tapply(rpsigdf$value, rpsigdf$replicate, mean)
                rep.nf <- min(avg_repsig) / avg_repsig

                ## apply normalization
                out <- mtgcovl
                for (spln in names(mtgcovl)) {
                    xxx = out[[spln]]
                    xxx$value <- xxx$value * rep.nf[desdf$replicate[desdf$sample == spln]]
                    out[[spln]] <- xxx
                }
                mtgcovl <- out
                remove(out)
            }

            return(mtgcovl)
        }

        ## for total RNA-seq
        {
            log_experiment <- desdf$RNA == 'total' & is.na(desdf$Sphase)
            mtgcovl <- cust_metagene_signal_fc(anngr, desdf=desdf[log_experiment,], upstream_length, downstream_length, extra_res, nbin, mtg_width, spi.seqdep.nf)

            ## plot the metagene of signal
            {
                ## reshape in dataframe
                gdata <- tibble(reshape2::melt(mtgcovl, id.vars=colnames(mtgcovl[[1]])))
                # remove(mtgcovl)

                colnames(gdata)[4] <- 'sample'
                gdata[c("condition", "replicate")] <- left_join(
                    gdata['sample'],
                    desdf,
                    by = 'sample'
                )[c("condition", "replicate")]
                gdata$replicate <- factor(gdata$replicate)

                fig <- ggplot(
                    gdata,
                    aes(
                        x = coord,
                        y = log10(value),
                        col = condition,
                        group = sample
                    )
                ) +
                    geom_line(aes(linetype = replicate), linewidth=2) + 
                    coord_cartesian(ylim=c(0, NA)) +
                    scale_color_manual(values=treat_colpal) +
                    geom_hline(yintercept=0, col='#000000') +
                    geom_vline(xintercept=0, col='#000000') +
                    geom_vline(xintercept=500, col='#000000') +
                    theme_txt
                
                figname <- file.path(figdir, 'metagene_signal_replHist_dwstrSL_total.pdf'); message(figname)
                pdf(figname, width = 2 * 7)
                print(fig)
                bouh <- dev.off()
            }
        }

        ## for 4sU RNA-seq
        {
            log_experiment <- desdf$RNA == 'nascent'
            mtgcovl <- cust_metagene_signal_fc(anngr, desdf=desdf[log_experiment,], upstream_length, downstream_length, extra_res, nbin, mtg_width, spi.seqdep.nf)

            ## plot the metagene of signal
            {
                ## reshape in dataframe
                gdata <- tibble(reshape2::melt(mtgcovl, id.vars=colnames(mtgcovl[[1]])))
                # remove(mtgcovl)

                colnames(gdata)[4] <- 'sample'
                gdata[c("condition", "replicate")] <- left_join(
                    gdata['sample'],
                    desdf,
                    by = 'sample'
                )[c("condition", "replicate")]
                gdata$replicate <- factor(gdata$replicate)

                fig <- ggplot(
                    gdata,
                    aes(
                        x = coord,
                        y = log10(value),
                        col = condition,
                        group = sample
                    )
                ) +
                    geom_line(aes(linetype = replicate), linewidth=2) + 
                    coord_cartesian(ylim=c(0, NA)) +
                    scale_color_manual(values=treat_colpal) +
                    geom_hline(yintercept=0, col='#000000') +
                    geom_vline(xintercept=0, col='#000000') +
                    geom_vline(xintercept=500, col='#000000') +
                    theme_txt
                
                figname <- file.path(figdir, 'metagene_signal_replHist_dwstrSL_nascent.pdf'); message(figname)
                pdf(figname, width = 2 * 7)
                print(fig)
                bouh <- dev.off()
            }
        }
    }

    ## clean memory
    bouh <- gc()


    #### get the foldchange per annotation
    {
        ## function to compute the metaplot signals
        cust_metagene_foldchange_fc <- function(anngr, desdf, upstream_length, downstream_length, extra_res, nbin, mtg_width, spi.seqdep.nf, test=FALSE)
        {
            tots <- list()
            covdfl <- list()
            if (test) {
                forvec <- c("D356T13", "D356T14", "D356T15", "D356T16") # TEMP
            } else {
                forvec <- desdf$sample
            }
    
            for (spln in forvec)
            {
                # spln <- "D356T17"

                ## load the coverage track (no 0 values)
                bwp <- file.path("Data", "ASF1_chaperone", "RNA-Seq", "coverage", paste0(spln, ".sorted.markdup.bw"))
                message(bwp)
                # file.exists(bwp)
                covgr <- rtracklayer::import.bw(bwp)

                ## recover signal along the metagene annotations
                {
                    ## find coverage blocks per replicative histone gene
                    ovl <- findOverlaps(covgr, anngr)
                    ovlt <- tapply(queryHits(ovl), subjectHits(ovl), identity)

                    ## compute the coverage track from the blocks
                    covl <- list()
                    tots[[spln]] <- list()
                    for (anni in 1:length(anngr)) {
                        # anni <- 1

                        ## extract the coverage blocks
                        scovgr <- covgr[ovlt[[anni]]]
                        start(scovgr)[1] <- start(anngr[anni])
                        end(scovgr)[length(scovgr)] <- end(anngr[anni])
                        scovgr$ordi <- as.numeric(1:length(scovgr))

                        ## identify the missing 0-value blocks
                        {
                            endv <- end(scovgr)
                            startv <- start(scovgr)
                            spacev <- tail(startv, n=-1) - head(endv, n=-1) - 1
                            log_zbl <- spacev != 0
                            if (any(log_zbl)) {
                                id_zbl <- which(log_zbl)

                                zgr <- scovgr[id_zbl]
                                zgr$score <- 0
                                zgr$ordi <- zgr$ordi + 0.5
                                end(zgr) <- start(scovgr)[id_zbl + 1] - 1
                                start(zgr) <- end(scovgr)[id_zbl] + 1

                                scovgr <- c(scovgr, zgr)
                                scovgr <- scovgr[order(scovgr$ordi)]
                            }
                        }


                        ## build coverage vector 
                        covv <- rep(scovgr$score, width(scovgr))
                        # print(c(tail(end(scovgr), n=1) - start(scovgr)[1] + 1, length(covv), width(anngr[anni])))
                        # if (tail(end(scovgr), n=1) - start(scovgr)[1] + 1 != length(covv)) {
                        #     message('longueurs différentes')
                        #     print(anni)
                        # }

                        ## normalize by spike-in controls and sequencing depth
                        covv <- spi.seqdep.nf$norm_factor[spi.seqdep.nf$sample == spln] * covv

                        ## get total signal in metagene normalized by gene length (before transformation in metagene)
                        tots[[spln]][[anngr$gene_name[anni]]] <- mean(covv[(upstream_length + 1):(upstream_length + width(anngr[anni]) - upstream_length - downstream_length)])

                        ## adjust by transcription sens
                        if (as.vector(strand(anngr[anni])) == '-') {
                            covv <- rev(covv)
                        }


                        ## transform to metagene scale
                        {
                            ## partition the signal
                            {
                                upv <- head(covv, upstream_length)
                                mtgv <- covv[(upstream_length + 1):(length(covv) - downstream_length)]
                                dwv <- tail(covv, downstream_length)
                                # all(c(upv, mtgv, dwv) == covv)
                            }

                            ## summarize the metagene part
                            {
                                binsz <- length(mtgv) / nbin
                                brds <- floor(binsz * 1:nbin)
                                ubrds <- c(1, head(brds, n=-1) + 1)

                                mtgbv <- apply(data.frame(u=ubrds, d=brds), 1, function(xxx, vec) {
                                    out <- vec[xxx[1]:xxx[2]]
                                    out <- mean(out)

                                    return(out)
                                }, vec=mtgv)
                            }
                        }

                        ## record the metagene track
                        covl[[anngr$gene_name[anni]]] <- c(upv, mtgbv, dwv)

                    }

                    ## reshape in data.frame
                    {
                        covdf <- as.data.frame(do.call(rbind, covl))

                        colnames(covdf) <- c(paste0('u', 1:upstream_length), paste0('m', 1:nbin), paste0('d', 1:downstream_length))
                        covdfl[[spln]] <- covdf
                    }
                }
            }


            ## normalize by total signal in the replicative histone gene in siCtrl condition
            {
                totsdf <- reshape2::melt(tots)
                colnames(totsdf) <- c('signal', 'gene_name', 'sample')
                # str(totsdf)
                stotsdf = totsdf
                totsdf$condition <- left_join(
                    totsdf['sample'],
                    desdf,
                    by = 'sample'
                )$condition
                totsdf <- pivot_wider(totsdf, names_from='gene_name', values_from='signal')

                ctots <- colSums(totsdf[totsdf$condition == "siCtrl", anngr$gene_name])
                ctots.log <- log10(ctots)
                ctots.log.avg <- mean(ctots.log)
                ctots.log.difftoavg <- ctots.log - ctots.log.avg
                sig.cds.nf <- 10^(-ctots.log.difftoavg)

                ncovdfl <- list()
                for (spln in names(covdfl)) {
                    covdf = covdfl[[spln]]

                    ncovdf <- as.data.frame(lapply(covdf, function(xxx, sig.cds.nf) {
                        out = xxx * sig.cds.nf
                        return(out)
                    }, sig.cds.nf))

                    ncovdfl[[spln]] <- ncovdf
                }
            }


            ## compute the fold change from siCtrl to siAsf1, by replicate
            {
                if (test) {
                    spl_by_rep_cond <- pivot_wider(desdf[desdf$RNA == 'total' & is.na(desdf$Sphase),], id_cols='replicate', names_from = 'condition', values_from = 'sample') # TEMP
                } else {
                    spl_by_rep_cond <- pivot_wider(desdf, id_cols='replicate', names_from = 'condition', values_from = 'sample')
                }

                ## compute ratio
                nfcl <- list()
                for (repn in spl_by_rep_cond$replicate)
                {
                    # repn <- 1
                    ncovdf_ctrl = ncovdfl[[spl_by_rep_cond$siCtrl[spl_by_rep_cond$replicate == repn]]]
                    ncovdf_treat = ncovdfl[[spl_by_rep_cond$siAsf1[spl_by_rep_cond$replicate == repn]]]

                    ## ratio treat / ctrl
                    nfcdf <- as.matrix(ncovdf_treat) / as.matrix(ncovdf_ctrl)

                    ## remove annotations NA or Inf values
                    log_nonnum <- apply(nfcdf, 1, function(xxx) {
                        any(is.na(xxx) | is.infinite(xxx))
                    })
                    nfcdf <- nfcdf[!log_nonnum,]

                    ## reshape
                    nfcdf <- as.data.frame(nfcdf)
                    nfcl[[repn]] <- nfcdf
                }
            }


            ## transform raw values in metagene track
            mtgfcl <- list()
            for (repn in names(nfcl))
            {
                nfcdf = nfcl[[repn]]

                ## compute meta fold change
                mtgfc <- unlist(lapply(nfcdf, mean))

                ## summarize extra regions in bins
                {
                    ## downstream region
                    {
                        downstream_length_fin <- floor(downstream_length / extra_res)
                        upstream_length_fin <- floor(upstream_length / extra_res)

                        dwm <- bin_signal(mtgfc[grepl('^d', names(nfcdf))], extra_res)
                        upm <- rev(bin_signal(rev(mtgfc[grepl('^u', names(nfcdf))]), extra_res))

                        ## adjust the names of the values
                        names(dwm) <- paste0('d', ((0:(downstream_length_fin - 1)) + 1) * extra_res)
                        names(upm) <- paste0('u', - ((0:(upstream_length_fin - 1)) + 1) * extra_res)

                        ## combine the regions
                        mtgfc <- c(upm, mtgfc[grepl('^m', names(mtgfc))], dwm)
                    }
                }

                ## reshape in dataframe
                mtgfc <- named_vec_to_df(mtgfc)

                ## set metagene coordinates
                {
                    coordv <- (-upstream_length_fin:-1) * extra_res
                    coordv <- c(coordv, 0:(nbin - 1) * (mtg_width / nbin))
                    coordv <- c(coordv, mtg_width + ((1:downstream_length_fin) - 1) * extra_res)
                    mtgfc$coord <- coordv
                }

                ## record metagene signal for the sample
                mtgfcl[[repn]] <- mtgfc
            }

            return(mtgfcl)
        }

        ## axes limits
        {
            ymin <- -0.7
            ymax <- 1.8
        }

        ## for total RNA-seq
        {
            log_experiment <- desdf$RNA == 'total' & is.na(desdf$Sphase)
            mtgfcl <- cust_metagene_foldchange_fc(anngr, desdf=desdf[log_experiment,], upstream_length, downstream_length, extra_res, nbin, mtg_width, spi.seqdep.nf)

            ## plot the metagene of signal
            {
                ## reshape in dataframe
                gdata <- tibble(reshape2::melt(mtgfcl, id.vars=colnames(mtgfcl[[1]])))
                # remove(mtgfcl)

                colnames(gdata)[4] <- 'replicate'
                gdata$replicate <- factor(gdata$replicate)
                
                fig <- ggplot(
                    gdata[!grepl('^u', gdata$name),],
                    aes(
                        x = coord,
                        y = log2(value),
                    )
                ) +
                    geom_line(aes(linetype = replicate), linewidth=2) + 
                    coord_cartesian(ylim=c(ymin, ymax)) +
                    # scale_color_manual(values=treat_colpal) +
                    geom_hline(yintercept=0, col='#000000') +
                    geom_vline(xintercept=0, col='#000000') +
                    geom_vline(xintercept=500, col='#000000') +
                    theme_txt
                
                figname <- file.path(figdir, 'metagene_foldChange_replHist_dwstrSL_total.pdf'); message(figname)
                pdf(figname, width = 1 * 7)
                print(fig)
                bouh <- dev.off()
            }
        }

        ## for 4sU RNA-seq
        {
            log_experiment <- desdf$RNA == 'nascent'
            mtgfcl <- cust_metagene_foldchange_fc(anngr, desdf=desdf[log_experiment,], upstream_length, downstream_length, extra_res, nbin, mtg_width, spi.seqdep.nf)

            ## plot the metagene of signal
            {
                ## reshape in dataframe
                gdata <- tibble(reshape2::melt(mtgfcl, id.vars=colnames(mtgfcl[[1]])))
                # remove(mtgfcl)

                colnames(gdata)[4] <- 'replicate'
                gdata$replicate <- factor(gdata$replicate)

                fig <- ggplot(
                    gdata[!grepl('^u', gdata$name),],
                    aes(
                        x = coord,
                        y = log2(value),
                    )
                ) +
                    geom_line(aes(linetype = replicate), linewidth=2) + 
                    coord_cartesian(ylim=c(ymin, ymax)) +
                    # scale_color_manual(values=treat_colpal) +
                    geom_hline(yintercept=0, col='#000000') +
                    geom_vline(xintercept=0, col='#000000') +
                    geom_vline(xintercept=500, col='#000000') +
                    theme_txt
                
                figname <- file.path(figdir, 'metagene_foldChange_replHist_dwstrSL_nascent.pdf'); message(figname)
                pdf(figname, width = 1 * 7)
                print(fig)
                bouh <- dev.off()
            }
        }
    }
}



####
