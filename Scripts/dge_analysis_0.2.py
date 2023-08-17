
import os
import re
import numpy
import patsy
import pandas
import seaborn as sns
import matplotlib.pyplot as plt

from collections import defaultdict
from matplotlib import colors, colorbar
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from sklearn.decomposition import PCA
from scipy.stats import pearsonr as pr
from scipy.cluster.hierarchy import fcluster, dendrogram

from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro

from utils.panels import MultiPanel

#differential expression analysis

edger = importr('edgeR')

pandas2ri.activate()


def qlf_test(r_fit,
             coef = -1,
             contrast = None):
    
    #F-test (by default on combined effect of all coefficients)
    if contrast is None:
        r_qlf = edger.glmQLFTest(r_fit, coef = coef)
    #F-test with contrasts
    else:
        r_qlf = edger.glmQLFTest(r_fit, contrast = contrast)

    #results sorted by FDR
    r_top = edger.topTags(r_qlf, 'all')
    
    return pandas2ri.ri2py(r_top[0]) 


def tmm_norm(counts,
             min_nonzero):
    
    #select genes with more than zero counts in at least n samples
    expr = counts[(counts > 0).sum(1) >= min_nonzero]
    
    #convert to R dataframe
    r_counts = pandas2ri.py2ri(expr)

    #DGEList object
    r_dgelist = edger.DGEList(counts = r_counts)
    
    #TMM-normalization
    r_dgelist = edger.calcNormFactors(r_dgelist, method = 'TMM')
    
    #normalized counts
    r_tmm = edger.cpm(r_dgelist)
    
    #convert to dataframe
    tmm = pandas2ri.ri2py(r_tmm)
    
    tmm = pandas.DataFrame(tmm, 
                           index = expr.index, 
                           columns = expr.columns)
    
    return tmm


def qlfit_edger(counts, 
                design, 
                formula,
                test_coef = -1,
                contrast = None,
                min_nonzero = 1):
    
    #design matrix
    dmatrix = patsy.dmatrix(formula, design)

    #raw counts
    counts = counts[design.index]
    
    #select genes with more than zero counts in at least n samples
    expr = counts[(counts > 0).sum(1) >= min_nonzero]
    
    #convert to R dataframe
    r_counts = pandas2ri.py2ri(expr)

    #DGEList object
    r_dgelist = edger.DGEList(counts = r_counts)
    
    #TMM-normalization
    r_dgelist = edger.calcNormFactors(r_dgelist, method = 'TMM')
    
    #estimate dispersion
    r_dgelist = edger.estimateDisp(r_dgelist, dmatrix)

    #fit GLM model
    r_fit = edger.glmQLFit(r_dgelist, dmatrix)
    
    #differential expression results
    results = qlf_test(r_fit, test_coef, contrast)
    
    #normalized counts
    r_tmm = edger.cpm(r_dgelist)
    
    #convert to dataframe
    tmm = pandas2ri.ri2py(r_tmm)
    
    tmm = pandas.DataFrame(tmm, index = expr.index, columns = expr.columns)

    results = results.join(tmm)

    return results, r_fit


################################################################################

#output path
if not os.path.exists('results/dge'):
    #create directory if necessary
    os.makedirs('results/dge')

#output path
if not os.path.exists('plots/dge'):
    #create directory if necessary
    os.makedirs('plots/dge/svg')

#output path
if not os.path.exists('results/gsea'):
    #create directory if necessary
    os.makedirs('results/gsea')

#global parameters for plotting

%matplotlib

sns.set(font = 'Liberation Sans', style = 'white')

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
                   
lblue, blue, dblue, ddblue, dddblue = '#92b7d1', '#3486bf', '#0a5c94', '#083b5e', '#0c1640'

lred, red, dred, ddred, dddred = '#cf8f8f', '#c46060', '#b02020', '#7a0c0c', '#4f0000'

################################################################################

#plot params

params = defaultdict(dict)

lblue, blue, dblue, ddblue, dddblue = '#92b7d1', '#3486bf', '#0a5c94', '#083b5e', '#0c1640'

lred, red, dred, ddred, dddred = '#cf8f8f', '#c46060', '#b02020', '#7a0c0c', '#4f0000'

params = defaultdict(dict)

for sample, (label, kd) in design_full[['label', 
                                     'knockdown']].iterrows():
    
    params['label'][sample] = label
    
    params['edgecolor'][sample] = red if kd else 'white'

    if 'nascent' not in label:
        
        params['color'][sample] = ddblue    
        
        boxcolor = ddred if kd else ddblue
        
    else:
        
        params['color'][sample] = dddblue

        boxcolor = dddred if kd else dddblue
            
    if '0.0h' in label:
    
        params['color'][sample] = lblue
        
        boxcolor = lred if kd else lblue  
        
    elif '2.0h' in label:
    
        params['color'][sample] = blue
        
        boxcolor = red if kd else blue 
        
    elif '5.0h' in label:
    
        params['color'][sample] = dblue
        
        boxcolor = dred if kd else dblue 
     
    params['boxcolor'][sample] = boxcolor 


################################################################################
################################################################################

#ERCC spike-in mix concentartions
conc = pandas.read_csv('data/annotation/cms_095046.txt', 
                       index_col = 'ERCC ID',
                       sep = '\t')

#gene and ERCC annotation
genes = pandas.read_csv('data/annotation/genes.GRCh38.95.csv',
                        index_col = 0)

spikes = pandas.read_csv('data/annotation/ERCC92.csv',
                         index_col = 0)

#gene name dictionaries
id2name = dict(genes['name'])

name2id = defaultdict(list)

for gene, name in id2name.items():
    
    name2id[name] += [gene] 

hist_genes = [gene for gene in id2name if id2name[gene].startswith('HIST')]

#experimental design (all samples)
design_full = pandas.read_csv('data/design.csv', index_col = 'sample')

design_full['group'] = design_full['label'].str[:-4]    

#all samples excluding asynchronous
sync = ~ design_full['time'].isnull()

#raw counts (including spike-in counts)
raw = pandas.read_csv('data/raw.csv.gz', index_col = 0)

#total number of mapped reads
tot = raw.sum()

#spike-in counts
raw_ercc = raw.loc[spikes.index]

#exclude spike-in counts
raw.drop(spikes.index, inplace = True)

#sample names dictionary
run2label = dict(design_full['label'])

#normalized counts

tmm_all = tmm_norm(raw[design_full.index], min_nonzero = 2)

logtmm_all = numpy.log2(1. + tmm_all)

relative_to_mean_all = logtmm_all.T - logtmm_all.T.mean()

################################################################################
#######################     DGE (async total)    ###############################
################################################################################

design_total = design_full[~ sync & (design_full['nascent'] == False)].copy()

design_total.sort_values('knockdown', inplace = True)

#normalized counts

tmm_total = tmm_norm(raw[design_total.index], min_nonzero = 2)

logtmm_total = numpy.log2(1. + tmm_total)

relative_to_mean_total = logtmm_total.T - logtmm_total.T.mean()

#differential expression analysis

thr = 0.05 #FDR 5%

results_total, r_fit_total = qlfit_edger(raw, 
                                         design_total, 
                                         formula = '~ knockdown',
                                         test_coef = ro.r('2'), 
                                         min_nonzero = 2)

results_total = results_total.join(genes)

sign_total = results_total[results_total['FDR'] < thr].index.tolist()

#save results

results_total.rename(columns = run2label).to_csv('results/dge/total.csv')

results_total.rename(columns = run2label).to_excel('results/dge/total.xlsx')

################################################################################
#                                 GSEA (total)                                 # 
################################################################################

#rank genes by effect i.e. logFC * - log10(PValue)

rnk_total = rank_genes(results_total)

#GO biological process (GSEA)

gsea_total = gsea_webgestalt(rnk_total, method = 'top', db = 'geneontology_Biological_Process_noRedundant') 

gsea_total.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_total.set_index('geneSet', inplace = True)

gsea_total.loc[:, 'names'] = rename(gsea_total['userId'])

gsea_total.to_csv('results/gsea/total_GO_BP.csv')

gsea_total.to_excel('results/gsea/total_GO_BP.xlsx')

#GO cellular component (GSEA)

gsea_cc_total = gsea_webgestalt(rnk_total, method = 'top', db = 'geneontology_Cellular_Component_noRedundant') 

gsea_cc_total.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_cc_total.set_index('geneSet', inplace = True)

gsea_cc_total.loc[:, 'names'] = rename(gsea_cc_total['userId'])

gsea_cc_total.to_csv('results/gsea/total_GO_CC.csv')

gsea_cc_total.to_excel('results/gsea/total_GO_CC.xlsx')

#KEGG (GSEA)

gsea_kegg_total = gsea_webgestalt(rnk_total, method = 'top', db = 'pathway_KEGG') 

gsea_kegg_total.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_kegg_total.set_index('geneSet', inplace = True)

gsea_kegg_total.loc[:, 'names'] = rename(gsea_kegg_total['userId'])

gsea_kegg_total.to_csv('results/gsea/total_KEGG.csv')

gsea_kegg_total.to_excel('results/gsea/total_KEGG.xlsx')

#WikiPathways (GSEA)

gsea_wiki_total = gsea_webgestalt(rnk_total, method = 'top', db = 'pathway_Wikipathway') 

gsea_wiki_total.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_wiki_total.set_index('geneSet', inplace = True)

gsea_wiki_total.loc[:, 'names'] = rename(gsea_wiki_total['userId'])

gsea_wiki_total.to_csv('results/gsea/total_WikiPathways.csv')

gsea_wiki_total.to_excel('results/gsea/total_WikiPathways.xlsx')

################################################################################
#######################     DGE (async nascent)    #############################
################################################################################

design_nascent = design_full[~ sync & (design_full['nascent'] == True)].copy()

design_nascent.sort_values('knockdown', inplace = True)

#normalized counts

tmm_nascent = tmm_norm(raw[design_nascent.index], min_nonzero = 2)

logtmm_nascent = numpy.log2(1. + tmm_nascent)

relative_to_mean_nascent = logtmm_nascent.T - logtmm_nascent.T.mean()

#differential expression analysis

thr = 0.05 #FDR 5%

results_nascent, r_fit_nascent = qlfit_edger(raw, 
                                             design_nascent, 
                                             formula = '~ knockdown',
                                             test_coef = ro.r('2'), 
                                             min_nonzero = 2)

results_nascent = results_nascent.join(genes)

sign_nascent = results_nascent[results_nascent['FDR'] < thr].index.tolist()

#save results

results_nascent.rename(columns = run2label).to_csv('results/dge/nascent.csv')

results_nascent.rename(columns = run2label).to_excel('results/dge/nascent.xlsx')

################################################################################
#                               GSEA (nascent)                                 # 
################################################################################

#rank genes by effect i.e. logFC * - log10(PValue)

rnk_nascent = rank_genes(results_nascent)

#GO biological process (GSEA)

gsea_nascent = gsea_webgestalt(rnk_nascent, method = 'top', db = 'geneontology_Biological_Process_noRedundant') 

gsea_nascent.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_nascent.set_index('geneSet', inplace = True)

gsea_nascent.loc[:, 'names'] = rename(gsea_nascent['userId'])

gsea_nascent.to_csv('results/gsea/nascent_GO_BP.csv')

gsea_nascent.to_excel('results/gsea/nascent_GO_BP.xlsx')

#GO cellular component (GSEA)

gsea_cc_nascent = gsea_webgestalt(rnk_nascent, method = 'top', db = 'geneontology_Cellular_Component_noRedundant') 

gsea_cc_nascent.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_cc_nascent.set_index('geneSet', inplace = True)

gsea_cc_nascent.loc[:, 'names'] = rename(gsea_cc_nascent['userId'])

gsea_cc_nascent.to_csv('results/gsea/nascent_GO_CC.csv')

gsea_cc_nascent.to_excel('results/gsea/nascent_GO_CC.xlsx')

#KEGG (GSEA)

gsea_kegg_nascent = gsea_webgestalt(rnk_nascent, method = 'top', db = 'pathway_KEGG') 

gsea_kegg_nascent.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_kegg_nascent.set_index('geneSet', inplace = True)

gsea_kegg_nascent.loc[:, 'names'] = rename(gsea_kegg_nascent['userId'])

gsea_kegg_nascent.to_csv('results/gsea/nascent_KEGG.csv')

gsea_kegg_nascent.to_excel('results/gsea/nascent_KEGG.xlsx')

#WikiPathways (GSEA)

gsea_wiki_nascent = gsea_webgestalt(rnk_nascent, method = 'top', db = 'pathway_Wikipathway') 

gsea_wiki_nascent.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_wiki_nascent.set_index('geneSet', inplace = True)

gsea_wiki_nascent.loc[:, 'names'] = rename(gsea_wiki_nascent['userId'])

gsea_wiki_nascent.to_csv('results/gsea/nascent_WikiPathways.csv')

gsea_wiki_nascent.to_excel('results/gsea/nascent_WikiPathways.xlsx')

################################################################################
#######################        DGE (sync GLM)    ###############################
################################################################################

#include synchronous samples only
design_sync = design_full[sync].copy()        

design_sync.sort_values(['knockdown', 'time'], inplace = True)

#normalized counts

tmm_sync = tmm_norm(raw[design_sync.index], min_nonzero = 2)

logtmm_sync = numpy.log2(1. + tmm_sync)

relative_to_mean_sync = logtmm_sync.T - logtmm_sync.T.mean()

#differential expression analysis

thr = 0.05 #FDR 5%

results_sync, r_fit_sync = qlfit_edger(raw, 
                                       design_sync, 
                                       formula = '~ C(time) * knockdown',
                                       #test combined effect of all coefficients 
                                       test_coef = ro.r('-1'), 
                                       min_nonzero = 2)

var_names = ['2.0h', '5.0h', 'siASF1', '2.0h:siASF1', '5.0h:siASF1']

results_sync.columns = var_names + list(results_sync.columns)[len(var_names):]

results_sync = results_sync.join(genes)

sign_sync = results_sync[results_sync['FDR'] < thr].index.tolist()

#sign = results[(results['FDR'] < thr) & (results[var_names].abs() > 1).any(1)].index.tolist()

#quasi-likelihood F-test (all coefficients)

diff = dict()

for i, name in enumerate(var_names):

    coef = i + 2
    
    diff[name] = qlf_test(r_fit_sync, coef = ro.r('%s' % coef))
    
    print(name, (diff[name]['FDR'] < thr).sum()) #FIXME

#save results

table_sync = results_sync[['logCPM', 'PValue', 'FDR']].join(genes)

for name in var_names:
    
    by_var = diff[name][['logFC', 'PValue', 'FDR']].loc[table_sync.index].copy()

    by_var.columns = '%s[logFC]' % name, '%s[PValue]' % name, '%s[FDR]' % name
    
    table_sync = table_sync.join(by_var)

table_sync.index.name = 'gene'

table_sync = table_sync.join(tmm_sync).rename(columns = run2label)

table_sync.to_csv('results/dge/sync.csv')

table_sync.to_excel('results/dge/sync.xlsx')

################################################################################
#                                GSEA (sync)                                   # 
################################################################################

#rank genes by effect i.e. logFC * - log10(PValue)

diff['siASF1']['PValue'].replace(0., 1e-200, inplace = True) #FIXME?

rnk_sync = rank_genes(diff['siASF1'])

#GO biological process (GSEA)

gsea_sync = gsea_webgestalt(rnk_sync, method = 'top', db = 'geneontology_Biological_Process_noRedundant') 

gsea_sync.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_sync.set_index('geneSet', inplace = True)

gsea_sync.loc[:, 'names'] = rename(gsea_sync['userId'])

gsea_sync.to_csv('results/gsea/sync_GO_BP.csv')

gsea_sync.to_excel('results/gsea/sync_GO_BP.xlsx')

#GO cellular component (GSEA)

gsea_cc_sync = gsea_webgestalt(rnk_sync, method = 'top', db = 'geneontology_Cellular_Component_noRedundant') 

gsea_cc_sync.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_cc_sync.set_index('geneSet', inplace = True)

gsea_cc_sync.loc[:, 'names'] = rename(gsea_cc_sync['userId'])

gsea_cc_sync.to_csv('results/gsea/sync_GO_CC.csv')

gsea_cc_sync.to_excel('results/gsea/sync_GO_CC.xlsx')

#KEGG (GSEA)

gsea_kegg_sync = gsea_webgestalt(rnk_sync, method = 'top', db = 'pathway_KEGG') 

gsea_kegg_sync.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_kegg_sync.set_index('geneSet', inplace = True)

gsea_kegg_sync.loc[:, 'names'] = rename(gsea_kegg_sync['userId'])

gsea_kegg_sync.to_csv('results/gsea/sync_KEGG.csv')

gsea_kegg_sync.to_excel('results/gsea/sync_KEGG.xlsx')

#WikiPathways (GSEA)

gsea_wiki_sync = gsea_webgestalt(rnk_sync, method = 'top', db = 'pathway_Wikipathway') 

gsea_wiki_sync.drop(['leadingEdgeId', 'plotPath'], 1, inplace = True)

gsea_wiki_sync.set_index('geneSet', inplace = True)

gsea_wiki_sync.loc[:, 'names'] = rename(gsea_wiki_sync['userId'])

gsea_wiki_sync.to_csv('results/gsea/sync_WikiPathways.csv')

gsea_wiki_sync.to_excel('results/gsea/sync_WikiPathways.xlsx')

################################################################################
#                                 misc                                         # 
################################################################################

#sets of differentially expressed genes

diff_set = dict()

for name in diff:
    
    diff_set[name] = diff[name][diff[name]['FDR'] < thr].index.tolist()
    
    diff_set[name] = set(diff_set[name])



plt.figure(); venn3([diff_set['siASF1'], diff_set['2.0h:siASF1'], diff_set['5.0h:siASF1']]) 

plt.figure(); venn3([diff_set['2.0h'], diff_set['5.0h'], diff_set['siASF1'] | diff_set['2.0h:siASF1'] | diff_set['5.0h:siASF1']])  

plt.figure(); venn3([diff_set['2.0h'], diff_set['5.0h'], diff_set['siASF1']])


