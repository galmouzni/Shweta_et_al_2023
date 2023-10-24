#R
{
    require(tidyverse)
    require(RColorBrewer)
    require(edgeR)
    require(rtracklayer)

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

        histgr <- rtracklayer::import.gff('./Data/histones.gff')
        hist_name_oldnew <- tibble(read.table('./Data/histone_geneNames_OldNew.tsv'))
        colnames(hist_name_oldnew) <- c('new', 'old')
        histgr_old <- histgr
        histgr_old$gene_name <- left_join(
            data.frame('new' = histgr$gene_name),
            hist_name_oldnew,
            by = 'new'
        )$old
    }


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

        ## adjust the table of read counts
        {
            ## combine the duplicated replicative histone genes
            for (dgi in 1:length(dupi.l)) {
                # dgi <- 1
                dgp <- dupi.l[[dgi]]
                
                stab <- atab[atab$gene_id %in% dgp,]
                # print(stab)

                sel_col <- !grepl('^(gene_id|PValue|FDR|name|description)', colnames(stab))
                stab[1, sel_col] <- stab[1, sel_col] + stab[2, sel_col]

                atab[atab$gene_id %in% dgp[1], sel_col] <- stab[1, sel_col]
                atab <- atab[!(atab$gene_id %in% dgp[2]),]
            }
        }
    }

    ## [Mix] make volcano plot siCtrl vs siAsf1 from asynchronized cells (total RNA-Seq) with histone genes (replicative and variants)
    {
        ## some computation
        {
            log_repl_hist <- grepl('^HIST', atab$name)
            log_vari_hist <- atab$gene_id %in% hctab$EnsemblGeneId[hctab$subClass == "replacement"]

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
        }

        ## print the plot
        {
            figname <- file.path(fig_dir, 'volcano_total.pdf')
            pdf(figname); message(figname)
            print(fig)
            bouh <- dev.off()
        }
    }
}




####
