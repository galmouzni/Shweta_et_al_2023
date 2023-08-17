#R
{
    require(tidyverse)
    require(RColorBrewer)
    require(edgeR)

    theme_txt <- theme(text = element_text(size = 20))
    old_theme <- theme_set(theme_bw())
}


## FUNCTIONS
{
    source('./Scripts/general/_functions_in_R.r')

    rgb_hex_to_fract <- function(cnn) {
        # cnn <- '#AE123A'
        cnn <- strsplit(sub('^#', '', cnn), '')[[1]]
        cnn <- c(paste0('0x', cnn[1], cnn[2]), paste0('0x', cnn[3], cnn[4]), paste0('0x', cnn[5], cnn[6]))
        cnn <- strtoi(cnn) / 256

        return(cnn)
    }

    color_grad_btw_hexRGB <- function (cone, ctwo, nnn) {
        # cone <- '#AE123A'
        # ctwo <- '#eeb85f'
        # nnn <- 63 # between 1 and 256

        cone <- rgb_hex_to_fract(cone)
        ctwo <- rgb_hex_to_fract(ctwo)

        nnn <- floor(nnn)
        cdiff <- ctwo - cone

        grad <- list(
            cone[1] + (0:(nnn - 1) / (nnn - 1)) * cdiff[1],
            cone[2] + (0:(nnn - 1) / (nnn - 1)) * cdiff[2],
            cone[3] + (0:(nnn - 1) / (nnn - 1)) * cdiff[3]
        )

        grad_vec <- rep("", length(grad[[1]]))
        for (iii in 1:length(grad[[1]])) {
            grad_vec[iii] <- rgb(grad[[1]][iii], grad[[2]][iii], grad[[3]][iii])
        }

        return(grad_vec)
    }

    color_gradient_multi <- function(col_lims, nnn=256) {
        # col_lims <- c(
        #     '#AE123A',
        #     '#eeb85f',
        #     '#f5ca63',
        #     '#e1bf67',
        #     '#2b5c8a'
        # )

        steps <- floor(1:(length(col_lims) - 1) * (nnn / (length(col_lims) - 1)))
        spaces <- steps - c(0, head(steps, n=-1))
        spaces[2:length(spaces)] <- spaces[2:length(spaces)] + 1

        final_pal <- list()
        for (iii in 1:(length(col_lims) - 1)) {
            final_pal <- c(final_pal, list(
                color_grad_btw_hexRGB(col_lims[iii], col_lims[iii + 1], spaces[iii])
            ))
        }

        if (length(col_lims) > 2) {
            for (iii in 2:length(final_pal)) {
                final_pal[[1]] <- c(final_pal[[1]], tail(final_pal[[iii]], n=-1))
            }
        }

        final_pal <- final_pal[[1]]

        return(final_pal)
    }

    color_palette_subset_by_values <- function(col_pal, vec, min_ext=min(vec, na.rm=TRUE), max_ext=max(vec, na.rm=TRUE)) {
        # col_pal <- rev(color_gradient_multi(col_lims=paste0('#', c("2e3193", "5073af", "86b2d1", "c2e0ec", "eff8e2", "fbeeae", "f5c17d", "e67f54", "c84133", "c84133")), nnn=256))
        # vec <- runif(10) * 2
        # min_ext = 0
        # max_ext = 2

        ## get the number of colors in the palette
        cpl <- length(col_pal)

        ## get min and max values of the subsetting vector
        min_val <- min(vec, na.rm=TRUE)
        max_val <- max(vec, na.rm=TRUE)

        ## find the relative positioning of the value min & max in the numerical color palette ranged
        pos_rel <- (c(min_val, max_val) - min_ext) / (max_ext - min_ext)

        ## subset the color_palette
        sub_col_pal <- col_pal[seq(from=floor(cpl * pos_rel[1]), to=floor(cpl * pos_rel[2]))]

        return(sub_col_pal)
    }
}


## MAIN
{
    ## general data
    {
        fig_dir <- "./Images/ASF1_chaperone/Shweta_reviews_230207"; dir.create(fig_dir, recursive=F, showWarnings=F)

        color_list <- list(
            'replicative' = '#A500CD', # '#CE00FF', # '#8d7493',
            'replicative_dark' = '#67007F',
            'replacement' = '#15D900', #'#11B300' # '#7bae6c'
            'replacement_dark' = '#3D7F26'
        )

        # col_pal_cust <- color_gradient_multi(col_lims=c('#AE123A', '#eeb85f', '#f5ca63', '#e1bf67', '#2b5c8a'), nnn=256)
        # col_pal_cust <- rev(color_gradient_multi(col_lims=paste0('#', c("003fa3","1e53a4","3b66a4","e9d8a6","ecba53","d47c2b","bb3e03","ab3015","9b2226")), nnn=256))
        # col_pal_cust <- rev(color_gradient_multi(col_lims=paste0('#', c("2e3193", "5073af", "86b2d1", "c2e0ec", "eff8e2", "fbeeae", "f5c17d", "e67f54", "c84133", "c84133")), nnn=256))
        # https://coolors.co/001219-005f73-0a9396-94d2bd-e9d8a6-ee9b00-ca6702-bb3e03-ae2012-9b2226
        col_pal_cust <- color_gradient_multi(col_lims=RColorBrewer::brewer.pal(n=11, name='RdYlBu'), nnn=256)
        
    }

    ## load list of genes
    {
        genep <- "./Data/GRCh38.104.sorted.geneOnly.gtf"
        genes <- rtracklayer::import.gff(genep)

        gene_conv_p <- "./Data/GRCh38_GRCh37_gene_names_correspondences.tsv"
        gene_conv <- read.delim(gene_conv_p, h=T, sep='\t') %>% tibble()

        gene_id_to_name <- as.list(
            gene_conv$gene_name.GRCh37
        )
        names(gene_id_to_name) <- genes$gene_id
    }

    ## load list of histone and chaperone
    {
        lhcp <- './Data/GtexTcgaGeneExpr/histone_chaperone_complete_list.csv'
        hctab <- read.delim(lhcp, h = T, sep = '\t', dec = ',') %>% tibble()
    }

    ## load the expression table from Alberto Gatto work
    {
        ## asynchronized cells
        {
            atabp <- "./Data/ASF1_chaperone/Alberto_Gatto_results/results/dge/total.tsv"
            atab <- read.delim(atabp, h = T, sep = '\t', dec = ',') %>% tibble()
            colnames(atab)[1] <- "gene_id"
            
            as_cnt_col <- c("control..1.", "control..2.", "siASF1..1.", "siASF1..2.")
            # atab[as_cnt_col]

            # atab_gn <- left_join(
            #     atab,
            #     gene_conv,
            #     by='gene_id'
            # )
            # atab_gn$log10_FDR <- -log10(atab_gn$FDR)
            # atab_gn <- atab_gn[c('gene_id', 'gene_name.GRCh38', 'gene_name.GRCh37', 'log10_FDR', colnames(atab_gn)[2:12])]
            # atab_gn[atab_gn$gene_name.GRCh38 == 'H3-3B' & !is.na(atab_gn$gene_name.GRCh38),]
        }

        ## G1/S synchronized cells
        if (FALSE) {
            # mtabp <- "../../Alberto Gatto-GA/Analyses/200707_new_rnaseq_results_shweta/analysis/results/dge/sync.xlsx"
            mtabp <- "./Data/ASF1_chaperone/Alberto_Gatto_results/results/dge/sync.tsv"
            # file.exists(mtabp)
            mtab <- read.delim(mtabp, h = T, sep = '\t', dec = ',') %>% tibble()
            colnames(mtab)[1] <- "gene_id"

            cnt_col <- c("control.0.0h..1.", "control.0.0h..2.", "control.2.0h..1.", "control.2.0h..2.", "control.5.0h..1.", "control.5.0h..2.", "siASF1.0.0h..1.", "siASF1.0.0h..2.", "siASF1.2.0h..1.", "siASF1.2.0h..2.", "siASF1.5.0h..1.", "siASF1.5.0h..2.")
            stat_col <- c("X2.0h.logFC.", "X2.0h.PValue.", "X2.0h.FDR.", "X5.0h.logFC.", "X5.0h.PValue.", "X5.0h.FDR.", "siASF1.logFC.", "siASF1.PValue.", "siASF1.FDR.", "X2.0h.siASF1.logFC.", "X2.0h.siASF1.PValue.", "X2.0h.siASF1.FDR.", "X5.0h.siASF1.logFC.", "X5.0h.siASF1.PValue.", "X5.0h.siASF1.FDR.")

            cltabp <- "./Data/ASF1_chaperone/Alberto_Gatto_results/results/clusters_6/gene_clusters_sync.tsv"
            cltab <- read.delim(cltabp, h = T, sep = '\t', dec = ',') %>% tibble()
            colnames(cltab)[1] <- "gene_id"
        }

        ## nascent RNAs
        if (FALSE) {
            rtabp <- "./Data/ASF1_chaperone/Alberto_Gatto_results/results/dge/nascent.tsv"
            rtab <- read.delim(rtabp, h = T, sep = '\t', dec = ',') %>% tibble()
            colnames(rtab)[1] <- "gene_id"

            rn_cnt_col <- c("control.nascent..1.", "control.nascent..2.", "siASF1.nascent..1.", "siASF1.nascent..2.")
            # rtab[rn_cnt_col]
        }
    }

    ## compute some general info
    {
        ## gene with significant differential expression in any comparison of synchronized cells
        {
            log_any_sig <- apply(mtab[stat_col[grepl('FDR', stat_col)]], 2, function(xxx) {
                xxx < 0.05
            })

            log_any_sig <- apply(log_any_sig, 1, any)
            # sum(log_any_sig)
            # 7706
        }
    }

    ## [Mix] make volcano plot siCtrl vs siAsf1 from asynchronized cells (total RNA-Seq) with histone genes (replicative and variants)
    {
        ## some computation
        {
            log_repl_hist <- grepl('^HIST', atab$name)
            log_vari_hist <- atab$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replacement"]

            ## check weird parts of the volcano plot (resolution to small ?)
            if (FALSE) {
                interc = 0.3
                slope = 1/5
                thresh_logfc = 5

                if (FALSE) {
                    fig <- ggplot(
                        atab,
                        aes(
                            x = logFC,
                            y = -log10(FDR)
                        )
                    ) +
                        geom_point(
                            col = '#CCCCCC'
                        ) +
                        geom_point(
                            data = atab[log_repl_hist,],
                            col = color_list$replicative,
                            size = hh_size,
                            alpha = hh_alpha
                        ) +
                        geom_point(
                            data = atab[log_vari_hist,],
                            col = color_list$replacement,
                            size = hh_size,
                            alpha = hh_alpha
                        ) +
                        geom_hline(yintercept = 1.3, linetype = 2) +
                        geom_vline(xintercept = 0, linetype = 2) +
                        geom_abline(intercept = interc, slope = slope, col = '#880000') +
                        geom_vline(xintercept = 5, col = '#880000') +
                        theme_txt

                    print(fig)
                }

                log_cust <- (-log10(atab$FDR) < atab$logFC * slope + interc) & (atab$logFC > thresh_logfc)
                quantile(atab[["control..1."]][log_cust], 1:9 / 10)
                quantile(atab[["control..2."]][log_cust], 1:9 / 10)

                ctrl_mean_cnt_cust <- rowMeans(atab[c("control..1.", "control..2.")][log_cust,])
                quantile(ctrl_mean_cnt_cust, 1:9 / 10)

                log_null_ctrl <- rowMeans(atab[c("control..1.", "control..2.")]) <= 0.1
                log_non_null_ctrl <- !log_null_ctrl

                log_cust <- (-log10(atab$FDR) > - atab$logFC * slope + interc) & (atab$logFC < - thresh_logfc / 2)
                quantile(atab[["siASF1..1."]][log_cust], 1:9 / 10)
                quantile(atab[["siASF1..2."]][log_cust], 1:9 / 10)

                log_null_siasf <- rowMeans(atab[c("siASF1..1.", "siASF1..2.")]) <= 0.1
                log_non_null_siasf <- !log_null_siasf

                log_non_null <- log_non_null_ctrl & log_non_null_siasf

                if (FALSE) {
                    fig <- ggplot(
                        atab[log_non_null,],
                        aes(
                            x = logFC,
                            y = -log10(FDR)
                        )
                    ) +
                        geom_point(
                            col = '#CCCCCC'
                        ) +
                        geom_point(
                            data = atab[log_repl_hist & log_non_null,],
                            col = color_list$replicative,
                            size = hh_size,
                            alpha = hh_alpha
                        ) +
                        geom_point(
                            data = atab[log_vari_hist & log_non_null,],
                            col = color_list$replacement,
                            size = hh_size,
                            alpha = hh_alpha
                        ) +
                        geom_hline(yintercept = 1.3, linetype = 2) +
                        geom_vline(xintercept = 0, linetype = 2) +
                        geom_abline(intercept = interc, slope = slope, col = '#880000') +
                        geom_vline(xintercept = 5, col = '#880000') +
                        theme_txt

                    print(fig)
                }
            }
            
        }

        ## [Plot] build the plot
        {
            hh_size = 5
            hh_alpha = 1 #0.80

            max_lfc <- max(atab$logFC) #[log_non_null])

            fig <- ggplot(
                atab, #[log_non_null,],
                aes(
                    x = logFC,
                    y = -log10(FDR)
                )
            ) +
                geom_point(
                    pch = 21,
                    fill = '#CCCCCC', # '#CCCCCC',
                    col = '#AAAAAA',
                    alpha = 1,
                    size = hh_size / 2
                ) +
                geom_point(
                    data = atab[log_repl_hist,], #[log_repl_hist & log_non_null,],
                    pch=21,
                    fill = color_list$replicative,
                    col = color_list$replicative_dark,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                geom_point(
                    data = atab[log_vari_hist,], #[log_vari_hist & log_non_null,],
                    pch=21,
                    fill = color_list$replacement,
                    col = color_list$replacement_dark,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                geom_hline(yintercept = 1.3, linetype = 2) +
                geom_vline(xintercept = 0, linetype = 2) +
                coord_cartesian(xlim=c(-max_lfc, max_lfc)) +
                theme_txt
            
            figname <- file.path(fig_dir, 'volcano_total.pdf')
            pdf(figname); message(figname)
            print(fig)
            bouh <- dev.off()
        }
    }

    ## [Mix] make volcano plot siCtrl vs siAsf1 from asynchronized cells (total RNA-Seq) with only the genes of group2 (see clusters based on the synchronized cells)
    if (FALSE) {
        ## some computation
        {
            log_group2 <- atab$gene_id %in% cltab$gene_id[cltab$cluster == 2]

            ## check weird parts of the volcano plot (resolution to small ?)
            {
                interc = 0.3
                slope = 1/5
                thresh_logfc = 5

                log_cust <- -log10(atab$FDR) < atab$logFC * slope + interc & atab$logFC > thresh_logfc
                quantile(atab[["control..1."]][log_cust], 1:9 / 10)
                quantile(atab[["control..2."]][log_cust], 1:9 / 10)

                ctrl_mean_cnt_cust <- rowMeans(atab[c("control..1.", "control..2.")][log_cust,])
                quantile(ctrl_mean_cnt_cust, 1:9 / 10)

                log_null_ctrl <- rowMeans(atab[c("control..1.", "control..2.")]) <= 0.1
                log_non_null_ctrl <- !log_null_ctrl

                log_cust <- -log10(atab$FDR) > - atab$logFC * slope + interc & atab$logFC < - thresh_logfc / 2
                quantile(atab[["siASF1..1."]][log_cust], 1:9 / 10)
                quantile(atab[["siASF1..2."]][log_cust], 1:9 / 10)

                log_null_siasf <- rowMeans(atab[c("siASF1..1.", "siASF1..2.")]) <= 0.1
                log_non_null_siasf <- !log_null_siasf

                log_non_null <- log_non_null_ctrl & log_non_null_siasf
            }
        }

        ## [Plot] build the plot
        {
            hh_size = 5
            hh_alpha = 0.80

            max_lfc <- max(atab$logFC[log_non_null])

            fig <- ggplot(
                atab[log_non_null,],
                aes(
                    x = logFC,
                    y = -log10(FDR)
                )
            ) +
                geom_point(
                    col = '#CCCCCC'
                ) +
                geom_point(
                    data = atab[log_group2 & log_non_null,],
                    col = color_list$replicative,
                    size = hh_size,
                    alpha = hh_alpha
                ) +
                geom_hline(yintercept = 1.3, linetype = 2) +
                geom_vline(xintercept = 0, linetype = 2) +
                coord_cartesian(xlim=c(-max_lfc, max_lfc)) +
                theme_txt
            
            figname <- file.path(fig_dir, 'volcano_total_group2HH.pdf')
            pdf(figname); message(figname)
            print(fig)
            bouh <- dev.off()
        }
    }

    ## [Mix] make the heatmap of expression ("z-scored") organized by clustering from synchronized cells
    if (FALSE) {
        ## prepare data
        {
            gdata <- left_join(
                cltab[c("gene_id", "cluster", "name")],
                mtab[c("gene_id", cnt_col)],
                by = "gene_id"
            )

            ## center and scale the counts
            # var_tot <- sd(unlist(gdata[cnt_col]))
            gdata[cnt_col] <- apply(gdata[cnt_col], 1, function(xxx) {
                (xxx - mean(xxx)) / (sd(xxx))
            }) %>% t() %>% as.data.frame()

            gdata <- gdata[order(gdata$cluster),]

            { ## order each gene cluster by hierarchical clustering on euclidian distance
                for (cii in unique(gdata$cluster)) {
                    # cii <- 2
                    log_cluster <- gdata$cluster == cii
                    sgdata <- as.matrix(gdata[log_cluster, cnt_col])
                    hcres <- hclust(as.dist(1 - cor(t(sgdata)) / 2), method='ward.D2')
                    # plot(hcres)
                    as.data.frame(gdata[log_cluster,])
                    gdata[log_cluster,] <- gdata[log_cluster,][hcres$order,]
                }
            }
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))

            ## check the cluster "2" obtained by Alberto only contains replicative histone
            # as.data.frame(gdata[gdata$cluster == 2,])

            gdata <- pivot_longer(gdata, cols=cnt_col, names_to = "condition", values_to = "count_scl_ctr")
            as.data.frame(gdata)[1:100,]
            cond_levs <- c(
                "control.0.0h..1.",
                "control.0.0h..2.",
                "control.2.0h..1.",
                "control.2.0h..2.",
                "control.5.0h..1.",
                "control.5.0h..2.",
                "siASF1.0.0h..1.",
                "siASF1.0.0h..2.",
                "siASF1.2.0h..1.",
                "siASF1.2.0h..2.",
                "siASF1.5.0h..1.",
                "siASF1.5.0h..2."
            )
            gdata$condition <- factor(gdata$condition, levels=rev(cond_levs))

            # gdata$count_scl_ctr <- (gdata$count - mean(gdata$count)) / sd(gdata$count)


            # gdata$count_scl_ctr <- gdata$count_scl_ctr * 20

            gdata$count_scl_ctr[gdata$count_scl_ctr > 2] <- 2
            gdata$count_scl_ctr[gdata$count_scl_ctr < -2] <- -2
            # summary(gdata$count_scl_ctr)

            # gdata$xxx <- 1/0
            # log_pos <- gdata$count_scl_ctr >= 0
            # log_neg <- !log_pos
            # gdata$xxx[log_pos] <- (gdata$count_scl_ctr[log_pos])^2
            # gdata$xxx[log_neg] <- -(-gdata$count_scl_ctr[log_neg])^2
            # # gdata$xxx <- gdata$xxx / (2^1.5)
            # gdata$count_scl_ctr <- gdata$xxx


            if (FALSE)
            {
                repart_plot(gdata$count_scl_ctr)

            }
        }

        ## [Plot] build plot
        {
            # col_pal_ext <- brewer.pal(n=3, name="RdYlBu")

            maxval <- max(abs(gdata$count_scl_ctr))
            col_xxx <- summary(floor((gdata$count_scl_ctr + maxval) / (2 * maxval) * (length(col_pal_cust) - 1) + 1))
            col_xxx <- col_xxx[1]:col_xxx[6]

            fig <- ggplot(
                gdata,
                aes(
                    x = gene_id,
                    y = condition,
                    fill = count_scl_ctr
                )
            ) +
                geom_tile() +
                scale_fill_gradientn(
                    # colors = rev(paletteer_c("grDevices::RdYlBu", 30))
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30))
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30))
                    colors = rev(col_pal_cust)[col_xxx],
                )
                # scale_fill_gradient2(
                #     low = col_pal_ext[3],
                #     mid = col_pal_ext[2],
                #     high = col_pal_ext[1],
                #     midpoint = 0
                # )

            figname <- file.path(fig_dir, "count_heatmap_synchronized_orderByCluster.pdf"); message(figname)
            pdf(figname, width = 4 * 7)
            print(fig)
            bouh <- dev.off()

            # figname <- file.path(fig_dir, "count_heatmap_synchronized_orderByCluster.png"); message(figname)
            # png(figname, width = 2 * 480, height = 0.5 * 480)
            # print(fig)
            # bouh <- dev.off()
        }

    }

    ## [Mix] make the heatmap of expression ("z-scored") organized by clustering from non-synchronized cells
    if (FALSE) {
        ## prepare data
        {
            gdata <- left_join(
                cltab[c("gene_id", "cluster", "name")],
                atab[c("gene_id", as_cnt_col)],
                by = "gene_id"
            )
            gdata <- gdata[order(gdata$cluster),]

            ## find NA genes
            log_na <- apply(gdata[, as_cnt_col], 1, function(xxx) {
                any(is.na(xxx))
            })
            log_nna <- !log_na

            ## scale and center
            gdata[as_cnt_col] <- apply(gdata[as_cnt_col], 1, function(xxx) {
                (xxx - mean(xxx)) / sd(xxx)
            }) %>% t() %>% as.data.frame()

            { ## order each gene cluster by hierarchical clustering on euclidian distance
                gdata$hcorder <- 0
                for (cii in unique(gdata$cluster)) {
                    log_cluster <- gdata$cluster == cii & log_nna
                    sgdata <- as.matrix(gdata[log_cluster, as_cnt_col])
                    hcres <- hclust(as.dist(1 - cor(t(sgdata)) / 2)) #, method='ward.D2')
                    # plot(hcres)
                    gdata$hcorder[log_cluster] <- hcres$order
                    gdata[log_cluster,] <- gdata[log_cluster,][gdata$hcorder[log_cluster],]
                }
            }
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))

            ## check the cluster "2" obtained by Alberto only contains replicative histone
            # as.data.frame(gdata[gdata$cluster == 2,])

            gdata <- pivot_longer(gdata, cols=as_cnt_col, names_to = "condition", values_to = "count_scl_ctr")
            cond_levs <- c(
                "control..1.",
                "control..2.",
                "siASF1..1.",
                "siASF1..2."
            )
            gdata$condition <- factor(gdata$condition, levels=rev(cond_levs))

            # gdata$count_scl_ctr <- (gdata$count - mean(gdata$count)) / sd(gdata$count)

            gdata$count_scl_ctr[gdata$count_scl_ctr > 2] <- 2
            gdata$count_scl_ctr[gdata$count_scl_ctr < -2] <- -2
        }

        ## [Plot] build plot
        {
            # col_pal_ext <- brewer.pal(n=3, name="RdYlBu")

            maxval <- max(abs(gdata$count_scl_ctr), na.rm=TRUE)
            col_xxx <- summary(floor((gdata$count_scl_ctr + maxval) / (2 * maxval) * (length(col_pal_cust) - 1) + 1))
            col_xxx <- col_xxx[1]:col_xxx[6]

            fig <- ggplot(
                gdata,
                aes(
                    x = gene_id,
                    y = condition,
                    fill = count_scl_ctr
                )
            ) +
                geom_tile() +
                scale_fill_gradientn(
                    # colors = rev(paletteer_c("grDevices::RdYlBu", 30))
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue Diverging", 30))
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30))
                    colors = rev(col_pal_cust)[col_xxx],
                )
                # scale_fill_gradient2(
                #     low = col_pal_ext[3],
                #     mid = col_pal_ext[2],
                #     high = col_pal_ext[1],
                #     midpoint = 0
                # )

            figname <- file.path(fig_dir, "count_heatmap_nonSynchro_orderByCluster.pdf"); message(figname)
            pdf(figname, width = 4 * 7)
            print(fig)
            bouh <- dev.off()
        }

    }

    ## [Mix] make the heatmap of expression ("z-scored") organized by clustering from nascent RNA signal
    if (FALSE) {
        ## prepare data
        {
            gdata <- left_join(
                cltab[c("gene_id", "cluster")],
                rtab[c("gene_id", rn_cnt_col)],
                by = "gene_id"
            )
            gdata <- gdata[order(gdata$cluster),]
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))

            gdata[rn_cnt_col] <- apply(gdata[rn_cnt_col], 1, function(xxx) {
                (xxx - mean(xxx)) / sd(xxx)
            }) %>% t() %>% as.data.frame()

            gdata <- pivot_longer(gdata, cols=rn_cnt_col, names_to = "condition", values_to = "count_scl_ctr")
            cond_levs <- c(
                "control.nascent..1.",
                "control.nascent..2.",
                "siASF1.nascent..1.",
                "siASF1.nascent..2."
            )
            gdata$condition <- factor(gdata$condition, levels=rev(cond_levs))

            # gdata$count_scl_ctr <- (gdata$count - mean(gdata$count)) / sd(gdata$count)

            gdata$count_scl_ctr[gdata$count_scl_ctr > 2] <- 2
            gdata$count_scl_ctr[gdata$count_scl_ctr < -2] <- -2
        }

        ## [Plot] build plot
        {
            maxval <- max(abs(gdata$count_scl_ctr), na.rm=T)
            col_xxx <- summary(floor((gdata$count_scl_ctr + maxval) / (2 * maxval) * (length(col_pal_cust) - 1) + 1))
            col_xxx <- col_xxx[1]:col_xxx[6]

            fig <- ggplot(
                gdata,
                aes(
                    x = gene_id,
                    y = condition,
                    fill = count_scl_ctr
                )
            ) +
                geom_tile() +
                scale_fill_gradientn(
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30))
                    colors = rev(col_pal_cust)[col_xxx]
                )

            figname <- file.path(fig_dir, "count_heatmap_nascent_orderByCluster.pdf")
            pdf(figname, width = 4 * 7); message(figname)
            print(fig)
            bouh <- dev.off()
        }
    }

    ## [Mix] make the heatmap of expression ("z-scored") for histone genes only (replicative, DEGs or not, variants) from nascent RNA signal
    if (FALSE) {
        ## prepare data
        {
            gdata <- left_join(
                cltab[c("gene_id", "cluster")],
                rtab[c("gene_id", rn_cnt_col)],
                by = "gene_id"
            )
            gdata <- gdata[order(gdata$cluster),]
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))

            gdata[rn_cnt_col] <- apply(gdata[rn_cnt_col], 1, function(xxx) {
                (xxx - mean(xxx)) / sd(xxx)
            }) %>% t() %>% as.data.frame()

            gdata <- gdata[gdata$gene_id %in% hctab$EnsemblGeneId[grepl('^Histone', hctab$Class)],]
            gdata$comment <- NA
            gdata$comment[gdata$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replicative"]] <- 1
            gdata$comment[gdata$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replacement"]] <- 3
            gdata$comment[gdata$cluster == 2] <- 2
            gdata <- gdata[order(gdata$comment),]
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))
            # gdata$name <- factor(gdata$name, levels=unique(gdata$name))
            # as.data.frame(gdata[1:5])

            gdata <- pivot_longer(gdata, cols=rn_cnt_col, names_to = "condition", values_to = "count_scl_ctr")
            cond_levs <- c(
                "control.nascent..1.",
                "control.nascent..2.",
                "siASF1.nascent..1.",
                "siASF1.nascent..2."
            )
            gdata$condition <- factor(gdata$condition, levels=rev(cond_levs))

            # gdata$count_scl_ctr <- (gdata$count - mean(gdata$count)) / sd(gdata$count)

            gdata$count_scl_ctr[gdata$count_scl_ctr > 2] <- 2
            gdata$count_scl_ctr[gdata$count_scl_ctr < -2] <- -2
        }

        ## [Plot] build plot
        {
            maxval <- max(abs(gdata$count_scl_ctr), na.rm=T)
            col_xxx <- summary(floor((gdata$count_scl_ctr + maxval) / (2 * maxval) * (length(col_pal_cust) - 1) + 1))
            col_xxx <- col_xxx[1]:col_xxx[6]

            fig <- ggplot(
                gdata,
                aes(
                    x = gene_id,
                    y = condition,
                    fill = count_scl_ctr
                )
            ) +
                geom_tile() +
                scale_fill_gradientn(
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30))
                    colors = rev(col_pal_cust)[col_xxx]
                )

            figname <- file.path(fig_dir, "count_heatmap_nascent_histonesOnly.pdf"); message(figname)
            pdf(figname, width = 1 * 7)
            print(fig)
            bouh <- dev.off()
        }
    }

    ## [Mix] make the heatmap of expression ("z-scored") for histone genes only (replicative, DEGs or not, variants) from synchronized cells
    if (FALSE) {
        ## prepare data
        {
            gdata <- left_join(
                cltab[c("gene_id", "name", "cluster")],
                mtab[c("gene_id", cnt_col)],
                by = "gene_id"
            )

            gdata[cnt_col] <- apply(gdata[cnt_col], 1, function(xxx) {
                (xxx - mean(xxx)) / sd(xxx)
            }) %>% t() %>% as.data.frame()

            gdata <- gdata[gdata$gene_id %in% hctab$EnsemblGeneId[grepl('^Histone', hctab$Class)],]
            gdata$comment <- NA
            gdata$comment[gdata$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replicative"]] <- 1
            gdata$comment[gdata$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replacement"]] <- 3
            gdata$comment[gdata$cluster == 2] <- 2
            gdata <- gdata[order(gdata$comment),]
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))
            gdata$name <- factor(gdata$name, levels=unique(gdata$name))
            # as.data.frame(gdata[1:5])

            gdata <- pivot_longer(gdata, cols=cnt_col, names_to = "condition", values_to = "count_scl_ctr")
            cond_levs <- c(
                "control.0.0h..1.",
                "control.0.0h..2.",
                "control.2.0h..1.",
                "control.2.0h..2.",
                "control.5.0h..1.",
                "control.5.0h..2.",
                "siASF1.0.0h..1.",
                "siASF1.0.0h..2.",
                "siASF1.2.0h..1.",
                "siASF1.2.0h..2.",
                "siASF1.5.0h..1.",
                "siASF1.5.0h..2."
            )
            gdata$condition <- factor(gdata$condition, levels=rev(cond_levs))

            # tot_sd <- sd(gdata$count_scl_ctr)
            # gdata$count_scl_ctr <- gdata$count_scl_ctr / tot_sd

            # gdata$count_scl_ctr <- (gdata$count - mean(gdata$count)) / sd(gdata$count)

            gdata$count_scl_ctr[gdata$count_scl_ctr > 2] <- 2
            gdata$count_scl_ctr[gdata$count_scl_ctr < -2] <- -2
        }

        ## [Plot] build plot
        {
            maxval <- max(abs(gdata$count_scl_ctr))
            col_xxx <- summary(floor((gdata$count_scl_ctr + maxval) / (2 * maxval) * (length(col_pal_cust) - 1) + 1))
            col_xxx <- col_xxx[1]:col_xxx[6]

            fig <- ggplot(
                gdata,
                aes(
                    x = name, #gene_id,
                    y = condition,
                    fill = count_scl_ctr
                )
            ) +
                geom_tile() +
                scale_fill_gradientn(
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30))
                    colors = rev(col_pal_cust)[col_xxx],
                    # values = c(-maxval, -maxval/2, 0, maxval/2, maxval)
                )
                #  +
                # theme(axis.text.x = element_text(angle = 45, hjust=1))

            figname <- file.path(fig_dir, "count_heatmap_synchronized_histonesOnly.pdf")
            pdf(figname, width = 7 * 1); message(figname)
            print(fig)
            bouh <- dev.off()
        }
    }

    ## [Mix] make the heatmap of expression ("z-scored") for histone genes only (replicative, DEGs or not, variants) from non-synchronized cells
    if (FALSE) {
        ## prepare data
        {
            gdata <- left_join(
                cltab[c("gene_id", "name", "cluster")],
                atab[c("gene_id", as_cnt_col)],
                by = "gene_id"
            )

            gdata[as_cnt_col] <- apply(gdata[as_cnt_col], 1, function(xxx) {
                (xxx - mean(xxx)) / sd(xxx)
            }) %>% t() %>% as.data.frame()

            gdata <- gdata[gdata$gene_id %in% hctab$EnsemblGeneId[grepl('^Histone', hctab$Class)],]
            gdata$comment <- NA
            gdata$comment[gdata$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replicative"]] <- 1
            gdata$comment[gdata$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replacement"]] <- 3
            gdata$comment[gdata$cluster == 2] <- 2
            gdata <- gdata[order(gdata$comment),]
            gdata$gene_id <- factor(gdata$gene_id, levels=unique(gdata$gene_id))
            gdata$name <- factor(gdata$name, levels=unique(gdata$name))
            # as.data.frame(gdata[1:5])

            gdata <- pivot_longer(gdata, cols=as_cnt_col, names_to = "condition", values_to = "count_scl_ctr")
            cond_levs <- c(
                "control..1.",
                "control..2.",
                "siASF1..1.",
                "siASF1..2."
            )
            gdata$condition <- factor(gdata$condition, levels=rev(cond_levs))

            # tot_sd <- sd(gdata$count_scl_ctr)
            # gdata$count_scl_ctr <- gdata$count_scl_ctr / tot_sd

            # gdata$count_scl_ctr <- (gdata$count - mean(gdata$count)) / sd(gdata$count)

            gdata$count_scl_ctr[gdata$count_scl_ctr > 2] <- 2
            gdata$count_scl_ctr[gdata$count_scl_ctr < -2] <- -2
        }

        ## [Plot] build plot
        {
            col_pal_cc <- color_palette_subset_by_values(col_pal_cust, gdata$count_scl_ctr, min_ext=-2, max_ext=2)

            fig <- ggplot(
                gdata,
                aes(
                    x = name, #gene_id,
                    y = condition,
                    fill = count_scl_ctr
                )
            ) +
                geom_tile() +
                scale_fill_gradientn(
                    # colors = rev(paletteer_c("ggthemes::Orange-Blue-White Diverging", 30))
                    colors = rev(col_pal_cc),
                    # values = c(-maxval, -maxval/2, 0, maxval/2, maxval)
                )
                #  +
                # theme(axis.text.x = element_text(angle = 45, hjust=1))

            figname <- file.path(fig_dir, "count_heatmap_nonSynchro_histonesOnly.pdf")
            pdf(figname, width = 7 * 1); message(figname)
            print(fig)
            bouh <- dev.off()
        }
    }
}




####
