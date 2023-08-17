
#functional enrichment analysis

webg = importr('WebGestaltR')


def ora_webgestalt(genes,
                   ref_genes, 
                   db = 'geneontology_Biological_Process_noRedundant',
                   id_type = 'ensembl_gene_id',
                   method = 'top',
                   min_size = 5):
    
    genes = list(genes)
    
    ref_genes = list(ref_genes)
    
    res = webg.WebGestaltR(enrichMethod = 'ORA',
                           organism = 'hsapiens',
                           enrichDatabase = db,
                           interestGene = genes,
                           interestGeneType = id_type,
                           referenceGene = ref_genes,
                           referenceGeneType = id_type,
                           sigMethod = method,
                           minNum = min_size,
                           isOutput = False,
                           nThread = 4)
    
    res = pandas2ri.ri2py(res)
    
    return res

    
################################################################################

#output path
if not os.path.exists('results/clusters_6'):
    #create directory if necessary
    os.makedirs('results/clusters_6')


#output path
if not os.path.exists('plots/clusters_6'):
    #create directory if necessary
    os.makedirs('plots/clusters_6/svg')


################################################################################
  
cg = sns.clustermap(relative_to_mean_sync[sign_sync], 
                    method = 'ward', 
                    cmap = 'RdYlBu_r',
                    figsize = (12.5, 2.5),
                    vmin = -2., vmax = 2.,
                    row_colors = [params['boxcolor'][run] for run in relative_to_mean_sync.index],
                    cbar_kws = {'ticks' : [-2., 0., 2.]})

col_dendro = cg.dendrogram_col


clusters = fcluster(col_dendro.linkage, t = 6, criterion = 'maxclust')

clusters = pandas.DataFrame(clusters,
                            index = sign,
                            columns = ['cluster'])

cluster_genes = dict(clusters.groupby('cluster').groups)


cg = sns.clustermap(relative_to_mean_sync[sign_sync], 
                    method = 'ward', 
                    cmap = 'RdYlBu_r',
                    figsize = (12.5, 2.5),
                    vmin = -2., vmax = 2.,
                    row_colors = [params['boxcolor'][run] for run in relative_to_mean_sync.index],
                    col_colors = [pal[i-1] for i in clusters['cluster']],
                    cbar_kws = {'ticks' : [-2., 0., 2.]})

################################################################################
################################################################################


clusters = clusters.sort_values(['cluster']).join(genes)

clusters.to_csv('results/clusters_6/gene_clusters_sync.csv')

clusters.to_excel('results/clusters_6/gene_clusters_sync.xlsx')


################################################################################
################################################################################

#ORA (Gene Ontology)

ora_res = []

for cl in cluster_genes:

    res = ora_webgestalt(cluster_genes[cl],
                         sign_sync,                    
                         db = ['geneontology_Biological_Process_noRedundant'], method = 'top')
    
    res.loc[:, 'cluster'] = cl
    
    #ora_plot(res)
    
    ora_res += [res]
                 
ora_res = pandas.concat(ora_res)

#ora_res = ora_res.set_index('geneSet')

ora_res.loc[:, 'names'] = rename(ora_res['userId'])

ora_res = ora_res[['cluster', 'geneSet', 'description', 'size', 'overlap', 'expect', 'enrichmentRatio', 'pValue', 'FDR', 'link', 'names']]

ora_res.to_csv('results/clusters_6/ora_GO_clusters_sync.csv', index = False)

ora_res.to_excel('results/clusters_6/ora_GO_clusters_sync.xlsx', index = False)

ora_res = ora_res.set_index('geneSet')


#ORA (KEGG)

ora_res_kegg = []

for cl in cluster_genes:

    res = ora_webgestalt(cluster_genes[cl],
                         sign_sync,                    
                         db = ['pathway_KEGG'], method = 'top')
    
    res.loc[:, 'cluster'] = cl
    
    #ora_plot(res)
    
    ora_res_kegg += [res]

ora_res_kegg = pandas.concat(ora_res_kegg)

#ora_res_kegg = ora_res.set_index('geneSet')


ora_res_kegg.loc[:, 'names'] = rename(ora_res_kegg['userId'])

ora_res_kegg = ora_res_kegg[['cluster', 'geneSet', 'description', 'size', 'overlap', 'expect', 'enrichmentRatio', 'pValue', 'FDR', 'link', 'names']]

ora_res_kegg.to_csv('results/clusters_6/ora_KEGG_clusters_sync.csv', index = False)

ora_res_kegg.to_excel('results/clusters_6/ora_KEGG_clusters_sync.xlsx', index = False)

ora_res_kegg = ora_res_kegg.set_index('geneSet')

################################################################################
################################################################################


def cluster_plot(subset, 
                 table,
                 cluster_name = '0', 
                 cluster_color = 'black',
                 maxc = 40,
                 filename = None, 
                 savefig = False):
    
    #upper panel (box plot)
    
    subset = [gene for gene in subset if gene in relative_to_mean_all]
    
    df = relative_to_mean_all[subset].join(design_full)

    df['time'] = [hours if hours >= 0. else (is_nascent and -2.) or -1. for hours, is_nascent in df[['time', 'nascent']].values]
    
    df = pandas.melt(df, id_vars = design_full.columns)
    
    #

    fig = MultiPanel(left = 0.9,
                     right = 0.3,
                     bottom = 0.3,
                     top = 0.35,
                     panel_height = 3.5, 
                     panel_width = 3.)

    ax = fig.axes[-1]

    sns.boxplot(x = 'time', 
                y = 'value', 
                hue = 'knockdown', 
                data = df, 
                ax = ax,
                notch = True,
                fliersize = 0.5,
                linewidth = 1.25,
                palette = [blue, red]) 

    ax.set_ylabel('Relative expression', fontsize = 14, fontweight = 'bold', labelpad = 15)

    plt.setp(ax.get_yticklabels(), fontsize = 12)

    ax.set_xlim(-0.75, 4.75)
    
    ax.set_ylim(-4.5, 4.5)

    ax.set_xlabel('')

    ax.set_xticklabels(['nascent', 'total', '0.0h', '2.0h', '5.0h'])

    plt.setp(ax.get_xticklabels(), fontsize = 12)

    plt.setp(ax.get_xticklabels()[2:], fontweight = 'bold')

    ax.legend().remove()        
    
    ax.axhline(0., color = 'black', linestyle = '--', linewidth = 1., alpha = 0.5, zorder = -10)    
    
    ax.spines['top'].set_visible(False)

    ax.spines['right'].set_visible(False)
    
    ax.set_title('Cluster %s' % cluster_name, fontsize = 14, fontweight = 'bold', color = cluster_color)
    
    #lower panel (ORA results)

    table = table.sort_values('FDR').head(5)
    
    table = table.sort_values('enrichmentRatio', ascending = True)
    
    fig.add_below(fig.panels[-1], height_frac = 1 /3.)

    ax = fig.panels[-1].axis    

    ax.set_xticks([])    
    ax.set_yticks([])
    
    [spine.set_visible(False) for loc, spine in ax.spines.items()]
    
    rows = table[['enrichmentRatio', 'FDR', 'description']].iterrows()

    for i, (geneset, (ratio, fdr, descr)) in enumerate(rows):
        
        sign = fdr < 0.05
        
        y = i 
        
        ax.text(0.0, y, 
                '%.1f' % ratio,
                ha = 'right',
                va = 'center', 
                fontweight = 'bold' if sign else 'medium',
                alpha = 1. if sign else 0.65,
                fontsize = 11, clip_on = False)

        ax.scatter(0.05, y + 0.05, 
                   s = 25,
                   facecolor = 'black' if sign else 'white',
                   edgecolor = 'black',
                   alpha = 1. if sign else 0.65)

        ax.text(0.1, y, 
                shorter_descr(descr, maxc = maxc),
                ha = 'left',
                va = 'center', 
                fontweight = 'medium',
                alpha = 1. if sign else 0.65,
                fontsize = 11, clip_on = False)
        
    ax.set_ylim(-1, 5)
    
    ax.set_xlim(0, 1)
    
    handles = [Line2D([0], [0], marker = 'o', lw = 0, mfc = 'black', mec = 'black')]

    labels = ['$\mathregular{FDR\/<}$5%']
    
    legend = ax.legend(handles,
                        labels,
                        title = '',
                        loc = 'lower right', 
                        bbox_to_anchor = (1.075, -0.2),# len(table) * 3 + 2),
                        bbox_transform = ax.transAxes,
                        handletextpad = -0.25, 
                        labelspacing = 0.1,
                        borderaxespad = 0.,
                        columnspacing = 0.25,
                        frameon = False,
                        markerscale = 0.8)
    
    plt.setp(legend.get_texts(), fontsize = 10.5, fontweight = 'bold')
    
    
    if savefig:
        
        if filename is None:
        
            filename = 'cluster_%s' % cluster_name
        
        fig.plot.savefig('plots/clusters_6/%s.png' % filename, dpi = 300)

        fig.plot.savefig('plots/clusters_6/svg/%s.svg' % filename)
        
    return




for cl in cluster_genes:
    
    fig = cluster_plot(cluster_genes[cl], 
                       ora_res[ora_res['cluster'] == cl],
                       cluster_name = str(cl),
                       cluster_color = pal[cl-1],
                       filename = 'cluster_%s' % cl,
                       savefig = True)



for cl in cluster_genes:
    
    fig = cluster_plot(cluster_genes[cl], 
                       ora_res_kegg[ora_res_kegg['cluster'] == cl],
                       cluster_name = str(cl),
                       cluster_color = pal[cl-1],
                       filename = 'cluster_%s_KEGG' % cl,
                       savefig = True)





