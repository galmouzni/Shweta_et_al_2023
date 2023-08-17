
import re
import seaborn as sns
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.cluster.hierarchy import fcluster, dendrogram
from utils.panels import MultiPanel
from matplotlib.patches import Rectangle
from matplotlib import colors, colorbar

if True:
    def add_cbar(ax,
                vmin,
                vmax,
                cmap = 'RdYlBu_r',
                fontsize = 10,
                fontweight = 'medium',
                orientation = 'vertical',
                ticks = [],
                ticklabels = [],
                tick_length = 1.5,
                label_pad = 20):
        #
        norm = colors.Normalize(vmin, vmax)
        #
        cbar = colorbar.ColorbarBase(ax=ax, 
                                    cmap=cmap,
                                    norm=norm,
                                    orientation = orientation)
        #
        ax.tick_params(pad = label_pad,
                    length = tick_length)
        #
        cbar.set_ticks(ticks)
        #
        if orientation == 'vertical':
            ax.set_yticklabels(ticklabels,
                            fontsize = fontsize,
                            fontweight = fontweight, 
                            ha = 'right')
        elif orientation == 'horizontal':
            ax.set_xticklabels(ticklabels,
                            fontsize = fontsize,
                            fontweight = fontweight, 
                            ha = 'center')
        #
        return ax, cbar


## general details
if True:
    ## colors and fonts
    lblue, blue, dblue, ddblue, dddblue = '#92b7d1', '#3486bf', '#0a5c94', '#083b5e', '#0c1640'
    lred, red, dred, ddred, dddred = '#cf8f8f', '#c46060', '#b02020', '#7a0c0c', '#4f0000'
    #
    #
    sns.set(font = 'Liberation Sans', style = 'white')
    #
    qualitative_pal = ('#ce6575', #red 
                    '#ae7baf', #purple
                    '#688db7', #blue
                    '#80b580', #green
                    '#e5cd5b', #yellow
                    '#82b3d6', #light blue
                    '#f5b46a', #orange
                    '#f48579', #bright red
                    '#b1d86f', #light green
                    '#7cc8cc', #green blue
                    '#9a93d6', #light purple
                    '#c4c2b6', #beige
                    '#a6a9bf') #silver
    #
    ## load list of genes
    genep = "./Data/GRCh38.104.sorted.geneOnly.gtf"
    genes = pd.read_csv(genep, sep='\t', skiprows=5, header=None, names=('ref', 'database', 'biotype', 'start', 'end', 'score', 'strand', 'other', 'attr'))
    genes['gene_id'] = [re.sub('".*$', '', re.sub('^.*gene_id\ "', '', sss)) for sss in genes.attr]
    genes['gene_name'] = [re.sub('".*$', '', re.sub('^.*gene_name\ "', '', sss)) for sss in genes.attr]
    #
    gene_conv_p = "./Data/GRCh38_GRCh37_gene_names_correspondences.tsv"
    gene_conv = pd.read_csv(gene_conv_p, sep='\t')
    gene_conv.index = gene_conv.gene_id
    gene_conv.rename(columns={'gene_name.GRCh37': 'gene_name_GRCh37', 'gene_name.GRCh38': 'gene_name_GRCh38'}, inplace=True)
    #
    ## load list of histone and chaperone
    lhcp = './Data/GtexTcgaGeneExpr/histone_chaperone_complete_list.csv'
    hctab = pd.read_csv(lhcp, sep = '\t', decimal = ',')
    hctab.index = hctab.EnsemblGeneId




##################################################################################
### Heatmap expression in synchro cells siASF1                                 ###
##################################################################################

## load the counts (already TMM normalized, but not transformed)
table_sync = pd.read_csv('./Data/ASF1_chaperone/Alberto_Gatto_results/results/dge/sync.tsv', sep = '\t', decimal = ',')

## prepare the count table
cnt_col = ['control 0.0h (1)', 'control 0.0h (2)', 'control 2.0h (1)', 'control 2.0h (2)', 'control 5.0h (1)', 'control 5.0h (2)', 'siASF1 0.0h (1)', 'siASF1 0.0h (2)', 'siASF1 2.0h (1)', 'siASF1 2.0h (2)', 'siASF1 5.0h (1)', 'siASF1 5.0h (2)']
tmm_sync = table_sync[cnt_col]
tmm_sync.index = table_sync.gene


## log-transform
logtmm_sync = np.log2(1. + tmm_sync)

## center the log counts
relative_to_mean_sync = logtmm_sync.T - logtmm_sync.T.mean()

## recover the significant genes (differentially expressed)
thr = 0.05 #FDR 5%
sign_sync = table_sync[table_sync['FDR'] < thr].gene.tolist()

## compute first clustering
cg = sns.clustermap(relative_to_mean_sync[sign_sync], 
                    method = 'ward', 
                    cmap = 'RdYlBu_r',
                    figsize = (12.5, 2.5),
                    vmin = -2., vmax = 2.,
                    # row_colors = [params['boxcolor'][run] for run in relative_to_mean_sync.index],
                    cbar_kws = {'ticks' : [-2., 0., 2.]})

## recover dendrogram and define clusters
col_dendro = cg.dendrogram_col
clusters = fcluster(col_dendro.linkage, t = 6, criterion = 'maxclust')
clusters = pd.DataFrame(clusters,
                        index = sign_sync,
                        columns = ['cluster'])
cluster_genes = dict(clusters.groupby('cluster').groups)

## [optional] recompute clustering with cluster information
cg = sns.clustermap(relative_to_mean_sync[sign_sync], 
                    method = 'ward', 
                    cmap = 'RdYlBu_r',
                    figsize = (12.5, 2.5),
                    vmin = -2., vmax = 2.,
                    # row_colors = [params['boxcolor'][run] for run in relative_to_mean_sync.index],
                    # col_colors = [pal[i-1] for i in clusters['cluster']],
                    cbar_kws = {'ticks' : [-2., 0., 2.]})

## reorder the indexes by the clustering
idx = cg.dendrogram_row.reordered_ind
clustered = cg.data2d
row_linkage = cg.dendrogram_row.linkage.copy()
idx_col = cg.dendrogram_col.reordered_ind


#params for plot
if True:
    params = defaultdict(dict)
    for i, sample in enumerate(cnt_col):
        # sample = cnt_col[0]
        #
        params['label'][sample] = re.sub('\ \([12]\)', '', sample)
        params['edgecolor'][sample] = red if re.search('ASF1', sample) else 'white'
        #
        if '0.0h' in sample:
            params['color'][sample] = lblue
            boxcolor = lred if re.search('ASF1', sample) else lblue  
        elif '2.0h' in sample:
            params['color'][sample] = blue
            boxcolor = red if re.search('ASF1', sample) else blue 
        elif '5.0h' in sample:
            params['color'][sample] = dblue
            boxcolor = dred if re.search('ASF1', sample) else dblue 
        #
        params['boxcolor'][sample] = boxcolor 


## build the heatmap of the significant genes by clusters
if True:
    ## create the scaffold
    fig = MultiPanel(panel_width = 0.9, 
                    panel_height = 2.5,
                    left = 0.1,
                    right = 0.5,
                    bottom = .25,
                    top = 0.25,
                    vspace = 0.05)
    #
    fig.add_after(fig.panels[-1], width_frac = 10. / 0.9)
    #
    ## add sample annotation boxes (on mid panel)
    ax = fig.panels[0].axis
    for i, run in enumerate(reversed(cnt_col)): #(design_sync.index[::-1]):
        box = Rectangle(xy = (0., i), 
                        width = 5, 
                        height = 1, 
                        facecolor = params['boxcolor'][run],
                        edgecolor = 'white',
                        linewidth = 1.1)
        ax.add_artist(box)
    #
    ax.set_ylim(0, len(cnt_col))
    ax.set_xlim(-0.1, 5.1)
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    for i, sample in enumerate(reversed(cnt_col)):
        label = params['label'][sample]
        ax.text(2.5, i + 0.45, label,
                ha = 'center',
                va = 'center',
                color = 'white',
                alpha = 1.,
                fontsize = 8.5, 
                fontweight = 'bold')
    #
    ## heatmap panel
    ax = fig.panels[1].axis
    ax.pcolormesh(clustered.loc[reversed(cnt_col)], cmap = 'RdYlBu_r', vmin = -2., vmax = 2., rasterized = True)
    ax.set_yticks([])
    ax.set_xticks([])
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    #
    ## colorbar for gradient legend
    fig.add_after(fig.panels[-1],
                height_frac = .5,
                width_frac = .04 * 2,
                space_frac = 2.5)
    #
    ax = fig.panels[-1].axis
    ax, cbar = add_cbar(ax=ax,
                        cmap = 'RdYlBu_r',
                        orientation = 'vertical',
                        vmin = -2.,
                        vmax = +2., 
                        ticks = [-2, -1, 0, +1, +2],
                        label_pad = 15,
                        ticklabels = ['$\mathregular{-2}$',
                                    '$\mathregular{-1}$',
                                    '$\mathregular{ 0}$',
                                    '$\mathregular{+1}$',
                                    '$\mathregular{+2}$'])
    #
    ## draw gene cluster annotations
    fig.add_below(fig.panels[1], width_frac = 10. / 0.9, height_frac = 0.05 / 2)
    ax = fig.panels[-1].axis
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    cluster_colors = clusters['cluster'][idx_col]
    cluster_colors = cluster_colors.values
    #
    cmap = colors.ListedColormap(qualitative_pal[:10])
    bounds = np.arange(1, 12, 1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    ax.pcolormesh(cluster_colors.reshape(-1, 1).T, cmap = cmap, vmin=1, vmax=10) #, norm = norm
    #
    prev_i = 0
    for i, c in enumerate(cluster_colors[:-1]):
        if c != cluster_colors[i+1]:
            #ax.axvline(i, color = 'black', lw = .75, alpha = .9)
            ax.axvline(i, color = 'white')
            ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
            #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
            prev_i = i
    #
    ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    #
    ax = fig.panels[1].axis
    for i, c in enumerate(cluster_colors[:-1]):
        if c != cluster_colors[i+1]:
            ax.axvline(i, color = 'black', lw = .75, alpha = .5)
            #ax.axvline(i, color = 'white')
    #
    ## write figure files
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_sync_all_no_row_cluster_6_clusters.png', dpi = 600)
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_sync_all_no_row_cluster_6_clusters.svg', dpi = 600)


## Focus on the histone variants
if True:
    ## focus on histone genes
    hc_rtms = relative_to_mean_sync.T.loc[[xxx for xxx in hctab.EnsemblGeneId if xxx in relative_to_mean_sync.columns]]
    #
    ## order by cluster
    hc_rtms = hc_rtms.join(clusters, how='left')
    hc_rtms = hc_rtms.join(hctab[['GeneOfficialSymbol', 'subClass', 'Class']])
    hc_rtms = hc_rtms.loc[[bool(re.search('Histone', xxx)) for xxx in list(hc_rtms.Class)]]
    #
    hc_rtms = hc_rtms.sort_values(['subClass', 'cluster'])
    hc_rtms = hc_rtms.join(gene_conv, how='left')
    #
    ## display ordered table of genes with clusters
    # xxx = pd.DataFrame({
    #     'cluster': list(hc_rtms.cluster),
    #     'gene_name_GRCh37': list(hc_rtms.gene_name_GRCh37),
    #     'gene_name_GRCh38': list(hc_rtms.gene_name_GRCh38)
    # })
    # print(xxx.to_string())

## build the heatmap for histone variants only
if True:
    ## create the scaffold
    panel_width = 0.4
    fig = MultiPanel(panel_width = 0.2, #0.9, 
                    panel_height = 2.5,
                    left = 0.1,
                    right = 0.5,
                    bottom = .25,
                    top = 0.25,
                    vspace = 0.05)
    #
    fig.add_after(fig.panels[-1], width_frac = 10. / 0.9)
    #
    ## add sample annotation boxes (on mid panel)
    ax = fig.panels[0].axis
    for i, run in enumerate(reversed(cnt_col)): #(design_sync.index[::-1]):
        box = Rectangle(xy = (0., i), 
                        width = 5, 
                        height = 1, 
                        facecolor = params['boxcolor'][run],
                        edgecolor = 'white',
                        linewidth = 1.1)
        ax.add_artist(box)
    #
    ax.set_ylim(0, len(cnt_col))
    ax.set_xlim(-0.1, 5.1)
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    for i, sample in enumerate(reversed(cnt_col)):
        label = params['label'][sample]
        ax.text(2.5, i + 0.45, label,
                ha = 'center',
                va = 'center',
                color = 'white',
                alpha = 1.,
                fontsize = 8.5, 
                fontweight = 'bold')
    #
    ## heatmap panel
    ax = fig.panels[1].axis
    ax.pcolormesh(pd.DataFrame(hc_rtms.T.loc[reversed(cnt_col)], dtype='float'), cmap = 'RdYlBu_r', vmin = -3., vmax = 3., rasterized = True)
    ax.set_yticks([])
    ax.set_xticks([])
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    #
    ## colorbar for gradient legend
    fig.add_after(fig.panels[-1],
                height_frac = .5,
                width_frac = .04 * 2,
                space_frac = 2.5)
    #
    ax = fig.panels[-1].axis
    ax, cbar = add_cbar(ax=ax,
                        cmap = 'RdYlBu_r',
                        orientation = 'vertical',
                        vmin = -3.,
                        vmax = +3., 
                        ticks = [-3, -2, -1, 0, +1, +2, +3],
                        label_pad = 15,
                        ticklabels = ['$\mathregular{-3}$',
                                    '$\mathregular{-2}$',
                                    '$\mathregular{-1}$',
                                    '$\mathregular{ 0}$',
                                    '$\mathregular{+1}$',
                                    '$\mathregular{+2}$',
                                    '$\mathregular{+3}$'])
    #
    # ## draw gene cluster annotations
    # fig.add_below(fig.panels[1], width_frac = 10. / 0.9, height_frac = 0.05 / 2)
    # ax = fig.panels[-1].axis
    # [spine.set_visible(False) for loc, spine in ax.spines.items()]
    # ax.set_yticks([])
    # ax.set_xticks([])
    # #
    # cluster_colors = clusters['cluster'][idx_col]
    # cluster_colors = cluster_colors.values
    # #
    # cmap = colors.ListedColormap(qualitative_pal[:10])
    # bounds = np.arange(1, 12, 1)
    # norm = colors.BoundaryNorm(bounds, cmap.N)
    # ax.pcolormesh(cluster_colors.reshape(-1, 1).T, cmap = cmap, vmin=1, vmax=10) #, norm = norm
    # #
    # prev_i = 0
    # for i, c in enumerate(cluster_colors[:-1]):
    #     if c != cluster_colors[i+1]:
    #         #ax.axvline(i, color = 'black', lw = .75, alpha = .9)
    #         ax.axvline(i, color = 'white')
    #         ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    #         #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    #         prev_i = i
    # #
    # ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    # #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    # #
    # ax = fig.panels[1].axis
    # for i, c in enumerate(cluster_colors[:-1]):
    #     if c != cluster_colors[i+1]:
    #         ax.axvline(i, color = 'black', lw = .75, alpha = .5)
    #         #ax.axvline(i, color = 'white')
    #
    ## write figure files
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_sync_all_no_row_histOnly.png', dpi = 600)
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_sync_all_no_row_histOnly.svg', dpi = 600)



##################################################################################
### Heatmap expression in non-synchro cells siASF1                             ###
##################################################################################

## load the counts (already TMM normalized, but not transformed)
table_async = pd.read_csv('./Data/ASF1_chaperone/Alberto_Gatto_results/results/dge/total.tsv', sep = '\t', decimal = ',')
table_async.rename(columns = {list(table_async)[0]:'gene'}, inplace=True)

## prepare the count table
as_cnt_col = ['control (1)', 'control (2)', 'siASF1 (1)', 'siASF1 (2)']
tmm_async = table_async[as_cnt_col]
tmm_async.index = table_async.gene

## harmonize with data from synchronized cells and clusters of genes
tmm_async = tmm_sync.join(tmm_async, how='left')[as_cnt_col]
tmm_async = tmm_async.loc[clustered.columns]


#plot params
if True:
    params = defaultdict(dict)
    for i, sample in enumerate(as_cnt_col):
        # sample = cnt_col[0]
        #
        params['label'][sample] = re.sub('\ \([12]\)', '', sample)
        params['edgecolor'][sample] = red if re.search('ASF1', sample) else 'white'
        #
        params['color'][sample] = blue
        params['boxcolor'][sample] = red if re.search('ASF1', sample) else lblue


## log-transform
logtmm_async = np.log2(1. + tmm_async)

## center the log counts
relative_to_mean_async = logtmm_async.T - logtmm_async.T.mean()

## build the heatmap
if True:
    ## create the scaffold
    fig = MultiPanel(panel_width = 0.9, 
                    panel_height = 0.8, #2.5,
                    left = 0.1,
                    right = 0.5,
                    bottom = .25,
                    top = 0.25,
                    vspace = 0.05)
    #
    fig.add_after(fig.panels[-1], width_frac = 10. / 0.9)
    #
    ## add sample annotation boxes (on mid panel)
    ax = fig.panels[0].axis
    for i, run in enumerate(reversed(as_cnt_col)): #(design_async.index[::-1]):
        box = Rectangle(xy = (0., i), 
                        width = 5, 
                        height = 1, 
                        facecolor = params['boxcolor'][run],
                        edgecolor = 'white',
                        linewidth = 1.1)
        ax.add_artist(box)
    #
    ax.set_ylim(0, len(as_cnt_col))
    ax.set_xlim(-0.1, 5.1)
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    for i, sample in enumerate(reversed(as_cnt_col)):
        label = params['label'][sample]
        ax.text(2.5, i + 0.45, label,
                ha = 'center',
                va = 'center',
                color = 'white',
                alpha = 1.,
                fontsize = 8.5, 
                fontweight = 'bold')
    #
    ## heatmap panel
    ax = fig.panels[1].axis
    ax.pcolormesh(relative_to_mean_async.loc[reversed(as_cnt_col)], cmap = 'RdYlBu_r', vmin = -2., vmax = 2., rasterized = True)
    ax.set_yticks([])
    ax.set_xticks([])
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    #
    ## colorbar for gradient legend
    fig.add_after(fig.panels[-1],
                height_frac = .5,
                width_frac = .04 * 2,
                space_frac = 2.5)
    #
    ax = fig.panels[-1].axis
    ax, cbar = add_cbar(ax=ax,
                        cmap = 'RdYlBu_r',
                        orientation = 'vertical',
                        vmin = -2.,
                        vmax = +2., 
                        ticks = [-2, -1, 0, +1, +2],
                        label_pad = 15,
                        ticklabels = ['$\mathregular{-2}$',
                                    '$\mathregular{-1}$',
                                    '$\mathregular{ 0}$',
                                    '$\mathregular{+1}$',
                                    '$\mathregular{+2}$'])
    #
    ## draw gene cluster annotations
    fig.add_below(fig.panels[1], width_frac = 10. / 0.9, height_frac = 0.05 / 2)
    ax = fig.panels[-1].axis
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    cluster_colors = clusters['cluster'][idx_col]
    cluster_colors = cluster_colors.values
    #
    cmap = colors.ListedColormap(qualitative_pal[:10])
    bounds = np.arange(1, 12, 1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    ax.pcolormesh(cluster_colors.reshape(-1, 1).T, cmap = cmap, vmin=1, vmax=10) #, norm = norm
    #
    prev_i = 0
    for i, c in enumerate(cluster_colors[:-1]):
        if c != cluster_colors[i+1]:
            #ax.axvline(i, color = 'black', lw = .75, alpha = .9)
            ax.axvline(i, color = 'white')
            ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
            #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
            prev_i = i
    #
    ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    #
    ax = fig.panels[1].axis
    for i, c in enumerate(cluster_colors[:-1]):
        if c != cluster_colors[i+1]:
            ax.axvline(i, color = 'black', lw = .75, alpha = .5)
            #ax.axvline(i, color = 'white')
    #
    ## write figure files
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_async_all_no_row_cluster_6_clusters.png', dpi = 600)
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_async_all_no_row_cluster_6_clusters.svg', dpi = 600)

## Focus on the histone variants
if True:
    ## focus on histone genes
    hc_rtmas = relative_to_mean_async.T.loc[[xxx for xxx in hctab.EnsemblGeneId if xxx in relative_to_mean_async.columns]]
    #
    ## order by cluster
    hc_rtmas = hc_rtmas.join(clusters, how='left')
    hc_rtmas = hc_rtmas.join(hctab[['GeneOfficialSymbol', 'subClass', 'Class']])
    # hc_rtmas = hc_rtmas.loc[[bool(re.search('Histone', xxx)) for xxx in list(hc_rtmas.Class)]]
    hc_rtmas = hc_rtms[['Class']].join(hc_rtmas, lsuffix='aaa')
    #
    hc_rtmas = hc_rtmas.sort_values(['subClass', 'cluster'])
    hc_rtmas = hc_rtmas.join(gene_conv, how='left')
    #
    ## display ordered table of genes with clusters
    # xxx = pd.DataFrame({
    #     'cluster': list(hc_rtmas.cluster),
    #     'gene_name_GRCh37': list(hc_rtmas.gene_name_GRCh37),
    #     'gene_name_GRCh38': list(hc_rtmas.gene_name_GRCh38)
    # })
    # print(xxx.to_string())

## build the heatmap for histone variants only
if True:
    ## create the scaffold
    panel_width = 0.4
    fig = MultiPanel(panel_width = 0.2, #0.9, 
                    panel_height = 0.8,
                    left = 0.1,
                    right = 0.5,
                    bottom = .25,
                    top = 0.25,
                    vspace = 0.05)
    #
    fig.add_after(fig.panels[-1], width_frac = 10. / 0.9)
    #
    ## add sample annotation boxes (on mid panel)
    ax = fig.panels[0].axis
    for i, run in enumerate(reversed(as_cnt_col)): #(design_sync.index[::-1]):
        box = Rectangle(xy = (0., i), 
                        width = 5, 
                        height = 1, 
                        facecolor = params['boxcolor'][run],
                        edgecolor = 'white',
                        linewidth = 1.1)
        ax.add_artist(box)
    #
    ax.set_ylim(0, len(as_cnt_col))
    ax.set_xlim(-0.1, 5.1)
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    for i, sample in enumerate(reversed(as_cnt_col)):
        label = params['label'][sample]
        ax.text(2.5, i + 0.45, label,
                ha = 'center',
                va = 'center',
                color = 'white',
                alpha = 1.,
                fontsize = 8.5, 
                fontweight = 'bold')
    #
    ## heatmap panel
    ax = fig.panels[1].axis
    ax.pcolormesh(pd.DataFrame(hc_rtmas.T.loc[reversed(as_cnt_col)], dtype='float'), cmap = 'RdYlBu_r', vmin = -2., vmax = 2., rasterized = True)
    ax.set_yticks([])
    ax.set_xticks([])
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    #
    ## colorbar for gradient legend
    fig.add_after(fig.panels[-1],
                height_frac = .5,
                width_frac = .04 * 2,
                space_frac = 2.5)
    #
    ax = fig.panels[-1].axis
    ax, cbar = add_cbar(ax=ax,
                        cmap = 'RdYlBu_r',
                        orientation = 'vertical',
                        vmin = -2.,
                        vmax = +2., 
                        ticks = [-2, -1, 0, +1, +2],
                        label_pad = 15,
                        ticklabels = ['$\mathregular{-2}$',
                                    '$\mathregular{-1}$',
                                    '$\mathregular{ 0}$',
                                    '$\mathregular{+1}$',
                                    '$\mathregular{+2}$'])
    #
    # ## draw gene cluster annotations
    # fig.add_below(fig.panels[1], width_frac = 10. / 0.9, height_frac = 0.05 / 2)
    # ax = fig.panels[-1].axis
    # [spine.set_visible(False) for loc, spine in ax.spines.items()]
    # ax.set_yticks([])
    # ax.set_xticks([])
    # #
    # cluster_colors = clusters['cluster'][idx_col]
    # cluster_colors = cluster_colors.values
    # #
    # cmap = colors.ListedColormap(qualitative_pal[:10])
    # bounds = np.arange(1, 12, 1)
    # norm = colors.BoundaryNorm(bounds, cmap.N)
    # ax.pcolormesh(cluster_colors.reshape(-1, 1).T, cmap = cmap, vmin=1, vmax=10) #, norm = norm
    # #
    # prev_i = 0
    # for i, c in enumerate(cluster_colors[:-1]):
    #     if c != cluster_colors[i+1]:
    #         #ax.axvline(i, color = 'black', lw = .75, alpha = .9)
    #         ax.axvline(i, color = 'white')
    #         ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    #         #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    #         prev_i = i
    # #
    # ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    # #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    # #
    # ax = fig.panels[1].axis
    # for i, c in enumerate(cluster_colors[:-1]):
    #     if c != cluster_colors[i+1]:
    #         ax.axvline(i, color = 'black', lw = .75, alpha = .5)
    #         #ax.axvline(i, color = 'white')
    #
    ## write figure files
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_async_all_no_row_histOnly.png', dpi = 600)
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_async_all_no_row_histOnly.svg', dpi = 600)





##################################################################################
### Heatmap Nascent RNA expression in non-synchro cells siASF1                 ###
##################################################################################

## load the counts (already TMM normalized, but not transformed)
table_nasync = pd.read_csv('./Data/ASF1_chaperone/Alberto_Gatto_results/results/dge/nascent.tsv', sep = '\t', decimal = ',')
table_nasync.rename(columns = {list(table_nasync)[0]:'gene'}, inplace=True)
table_nasync.rename(columns = {'control nascent (1)':'control (1)', 'control nascent (2)':'control (2)', 'siASF1 nascent (1)':'siASF1 (1)', 'siASF1 nascent (2)':'siASF1 (2)'}, inplace=True)

## prepare the count table
nas_cnt_col = ['control (1)', 'control (2)', 'siASF1 (1)', 'siASF1 (2)']
tmm_nasync = table_nasync[nas_cnt_col]
tmm_nasync.index = table_nasync.gene

## harmonize with data from synchronized cells and clusters of genes
tmm_nasync = tmm_sync.join(tmm_nasync, how='left')[nas_cnt_col]
tmm_nasync = tmm_nasync.loc[clustered.columns]


#plot params
if True:
    params = defaultdict(dict)
    for i, sample in enumerate(nas_cnt_col):
        # sample = cnt_col[0]
        #
        params['label'][sample] = re.sub('\ \([12]\)', '', sample)
        params['edgecolor'][sample] = red if re.search('ASF1', sample) else 'white'
        #
        params['color'][sample] = blue
        params['boxcolor'][sample] = red if re.search('ASF1', sample) else lblue


## log-transform
logtmm_nasync = np.log2(1. + tmm_nasync)

## center the log counts
relative_to_mean_nasync = logtmm_nasync.T - logtmm_nasync.T.mean()

## build the heatmap
if True:
    ## create the scaffold
    fig = MultiPanel(panel_width = 0.9, 
                    panel_height = 0.8, #2.5,
                    left = 0.1,
                    right = 0.5,
                    bottom = .25,
                    top = 0.25,
                    vspace = 0.05)
    #
    fig.add_after(fig.panels[-1], width_frac = 10. / 0.9)
    #
    ## add sample annotation boxes (on mid panel)
    ax = fig.panels[0].axis
    for i, run in enumerate(reversed(nas_cnt_col)): #(design_nasync.index[::-1]):
        box = Rectangle(xy = (0., i), 
                        width = 5, 
                        height = 1, 
                        facecolor = params['boxcolor'][run],
                        edgecolor = 'white',
                        linewidth = 1.1)
        ax.add_artist(box)
    #
    ax.set_ylim(0, len(nas_cnt_col))
    ax.set_xlim(-0.1, 5.1)
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    for i, sample in enumerate(reversed(nas_cnt_col)):
        label = params['label'][sample]
        ax.text(2.5, i + 0.45, label,
                ha = 'center',
                va = 'center',
                color = 'white',
                alpha = 1.,
                fontsize = 8.5, 
                fontweight = 'bold')
    #
    ## heatmap panel
    ax = fig.panels[1].axis
    ax.pcolormesh(relative_to_mean_nasync.loc[reversed(nas_cnt_col)], cmap = 'RdYlBu_r', vmin = -2., vmax = 2., rasterized = True)
    ax.set_yticks([])
    ax.set_xticks([])
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    #
    ## colorbar for gradient legend
    fig.add_after(fig.panels[-1],
                height_frac = .5,
                width_frac = .04 * 2,
                space_frac = 2.5)
    #
    ax = fig.panels[-1].axis
    ax, cbar = add_cbar(ax=ax,
                        cmap = 'RdYlBu_r',
                        orientation = 'vertical',
                        vmin = -2.,
                        vmax = +2., 
                        ticks = [-2, -1, 0, +1, +2],
                        label_pad = 15,
                        ticklabels = ['$\mathregular{-2}$',
                                    '$\mathregular{-1}$',
                                    '$\mathregular{ 0}$',
                                    '$\mathregular{+1}$',
                                    '$\mathregular{+2}$'])
    #
    ## draw gene cluster annotations
    fig.add_below(fig.panels[1], width_frac = 10. / 0.9, height_frac = 0.05 / 2)
    ax = fig.panels[-1].axis
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    cluster_colors = clusters['cluster'][idx_col]
    cluster_colors = cluster_colors.values
    #
    cmap = colors.ListedColormap(qualitative_pal[:10])
    bounds = np.arange(1, 12, 1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    ax.pcolormesh(cluster_colors.reshape(-1, 1).T, cmap = cmap, vmin=1, vmax=10) #, norm = norm
    #
    prev_i = 0
    for i, c in enumerate(cluster_colors[:-1]):
        if c != cluster_colors[i+1]:
            #ax.axvline(i, color = 'black', lw = .75, alpha = .9)
            ax.axvline(i, color = 'white')
            ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
            #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
            prev_i = i
    #
    ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    #
    ax = fig.panels[1].axis
    for i, c in enumerate(cluster_colors[:-1]):
        if c != cluster_colors[i+1]:
            ax.axvline(i, color = 'black', lw = .75, alpha = .5)
            #ax.axvline(i, color = 'white')
    #
    ## write figure files
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_nasync_all_no_row_cluster_6_clusters.png', dpi = 600)
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_nasync_all_no_row_cluster_6_clusters.svg', dpi = 600)

## Focus on the histone variants
if True:
    ## focus on histone genes
    hc_rtmnas = relative_to_mean_nasync.T.loc[[xxx for xxx in hctab.EnsemblGeneId if xxx in relative_to_mean_nasync.columns]]
    #
    ## order by cluster
    hc_rtmnas = hc_rtmnas.join(clusters, how='left')
    hc_rtmnas = hc_rtmnas.join(hctab[['GeneOfficialSymbol', 'subClass', 'Class']])
    hc_rtmnas = hc_rtmnas.loc[[bool(re.search('Histone', xxx)) for xxx in list(hc_rtmnas.Class)]]
    #
    hc_rtmnas = hc_rtmnas.sort_values(['subClass', 'cluster'])
    hc_rtmnas = hc_rtmnas.join(gene_conv, how='left')
    #
    ## display ordered table of genes with clusters
    # xxx = pd.DataFrame({
    #     'cluster': list(hc_rtmnas.cluster),
    #     'gene_name_GRCh37': list(hc_rtmnas.gene_name_GRCh37),
    #     'gene_name_GRCh38': list(hc_rtmnas.gene_name_GRCh38)
    # })
    # print(xxx.to_string())

## build the heatmap for histone variants only
if True:
    ## create the scaffold
    panel_width = 0.4
    fig = MultiPanel(panel_width = 0.2, #0.9, 
                    panel_height = 0.8,
                    left = 0.1,
                    right = 0.5,
                    bottom = .25,
                    top = 0.25,
                    vspace = 0.05)
    #
    fig.add_after(fig.panels[-1], width_frac = 10. / 0.9)
    #
    ## add sample annotation boxes (on mid panel)
    ax = fig.panels[0].axis
    for i, run in enumerate(reversed(as_cnt_col)): #(design_sync.index[::-1]):
        box = Rectangle(xy = (0., i), 
                        width = 5, 
                        height = 1, 
                        facecolor = params['boxcolor'][run],
                        edgecolor = 'white',
                        linewidth = 1.1)
        ax.add_artist(box)
    #
    ax.set_ylim(0, len(as_cnt_col))
    ax.set_xlim(-0.1, 5.1)
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    ax.set_yticks([])
    ax.set_xticks([])
    #
    for i, sample in enumerate(reversed(as_cnt_col)):
        label = params['label'][sample]
        ax.text(2.5, i + 0.45, label,
                ha = 'center',
                va = 'center',
                color = 'white',
                alpha = 1.,
                fontsize = 8.5, 
                fontweight = 'bold')
    #
    ## heatmap panel
    ax = fig.panels[1].axis
    ax.pcolormesh(pd.DataFrame(hc_rtmnas.T.loc[reversed(as_cnt_col)], dtype='float'), cmap = 'RdYlBu_r', vmin = -2., vmax = 2., rasterized = True)
    xxx = hc_rtmnas.T.loc[reversed(as_cnt_col)]
    ax.set_yticks([])
    ax.set_xticks([])
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    #
    ## colorbar for gradient legend
    fig.add_after(fig.panels[-1],
                height_frac = .5,
                width_frac = .04 * 2,
                space_frac = 2.5)
    #
    ax = fig.panels[-1].axis
    ax, cbar = add_cbar(ax=ax,
                        cmap = 'RdYlBu_r',
                        orientation = 'vertical',
                        vmin = -2.,
                        vmax = +2., 
                        ticks = [-2, -1, 0, +1, +2],
                        label_pad = 15,
                        ticklabels = ['$\mathregular{-2}$',
                                    '$\mathregular{-1}$',
                                    '$\mathregular{ 0}$',
                                    '$\mathregular{+1}$',
                                    '$\mathregular{+2}$'])
    #
    # ## draw gene cluster annotations
    # fig.add_below(fig.panels[1], width_frac = 10. / 0.9, height_frac = 0.05 / 2)
    # ax = fig.panels[-1].axis
    # [spine.set_visible(False) for loc, spine in ax.spines.items()]
    # ax.set_yticks([])
    # ax.set_xticks([])
    # #
    # cluster_colors = clusters['cluster'][idx_col]
    # cluster_colors = cluster_colors.values
    # #
    # cmap = colors.ListedColormap(qualitative_pal[:10])
    # bounds = np.arange(1, 12, 1)
    # norm = colors.BoundaryNorm(bounds, cmap.N)
    # ax.pcolormesh(cluster_colors.reshape(-1, 1).T, cmap = cmap, vmin=1, vmax=10) #, norm = norm
    # #
    # prev_i = 0
    # for i, c in enumerate(cluster_colors[:-1]):
    #     if c != cluster_colors[i+1]:
    #         #ax.axvline(i, color = 'black', lw = .75, alpha = .9)
    #         ax.axvline(i, color = 'white')
    #         ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    #         #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    #         prev_i = i
    # #
    # ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'medium', color = 'black', ha = 'center', va = 'center', clip_on = False)
    # #ax.text(prev_i + (i - prev_i) / 2., -2, str(c), fontsize = 11, fontweight = 'bold', color = pal[c-1], ha = 'center', va = 'center', clip_on = False)
    # #
    # ax = fig.panels[1].axis
    # for i, c in enumerate(cluster_colors[:-1]):
    #     if c != cluster_colors[i+1]:
    #         ax.axvline(i, color = 'black', lw = .75, alpha = .5)
    #         #ax.axvline(i, color = 'white')
    #
    ## write figure files
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_nasync_all_no_row_histOnly.png', dpi = 600)
    fig.plot.savefig('./Images/ASF1_chaperone/Shweta_reviews_230207/4_heatmap_nasync_all_no_row_histOnly.svg', dpi = 600)



