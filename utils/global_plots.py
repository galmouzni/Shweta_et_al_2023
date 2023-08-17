#!/usr/bin/env python

import re
import numpy
import pandas
#import seaborn
#from matplotlib import pyplot
from scipy.stats import pearsonr, spearmanr

import seaborn as sns
import matplotlib.pyplot as plt

#global parameters 

sns.set(font = 'Liberation Sans', style = 'white')

blue = '#5b87b2ff'

purple = '#835d86ff'

white_purple = '#d1ccd1'

light_purple = '#b59db5'

dark_purple = '#664868ff'

green = '#559f53ff'

white_green = '#d3e2c9'

light_green = '#c5dbc5'

dark_green = '#1c6b23'

yellow = '#eac462'

greyish = '#e3e3e3ff'

blueish = '#5e7a91'

red = '#c44e52'

darker_blue = '#4c72b0'

purple_cmap = sns.blend_palette([white_purple, 
                                 light_purple, 
                                 purple, 
                                 dark_purple], as_cmap = True)

green_cmap = sns.blend_palette([white_green, 
                                light_green, 
                                green, 
                                dark_green], as_cmap = True)

#functions

def add_mask(ax,
             alpha = 0.5,
             linewidth = 1.1,
             linestyle = '-',
             line_color = greyish,
             mask_color = 'white'):
    
    #x-axis range
    xlim = ax.get_xlim()
    
    #y-axis lowerbound
    ymin = ax.get_ylim()[0]
    
    #maximum z-order in subplot
    zmax = max([child.zorder for child in ax.get_children()])
    
    #add horizontal line
    ax.axhline(linestyle = linestyle,
               linewidth = linewidth,
               color = line_color,
               zorder = zmax)
    
    #mask negative values (from ymin to zero)
    ax.fill_between(x = xlim, 
                    y1 = ymin, 
                    y2 = 0,
                    linewidth = 0.,
                    color = mask_color,
                    zorder = zmax,
                    alpha = alpha)
    
    for loc in ax.spines:
        #update z-order of spines
        ax.spines[loc].set_zorder(zmax + 1)
    
    return ax  


def fill_cov(ax, 
             cov,
             fill_color,
             line_color,
             zorder = 0,
             alpha = 1.,
             linewidth = 0.5):
    
    #x-axis range
    length = len(cov)
    
    #x-axis values (bin positions)
    pos = numpy.arange(length)
    
    #plot line
    ax.plot(x = pos, 
            y = cov, 
            color = line_color, 
            linewidth = linewidth,
            zorder = zorder - 1,
            alpha = alpha)
    
    #fill area
    ax.fill_between(x = pos, 
                    y1 = cov, 
                    color = fill_color,
                    linewidth = linewidth,
                    zorder = zorder - 1,
                    alpha = alpha)

    #set x-axis bounds
    ax.set_xlim(0, length)

    #hide x-axis ticklabels
    ax.set_xticklabels([])

    for loc in ['top', 'right', 'bottom']:
        #hide spines
        ax.spines[loc].set_visible(False)
    
    return ax


def fill_domains(ax, 
                 domains,
                 chromosome,
                 start,
                 end,
                 color,
                 ymin = 1.0,
                 ymax = 1.1,
                 alpha = 0.9,
                 binsize = 10000):

    bin_start = start // binsize

    bin_end = end // binsize + 1
    
    length = bin_end - bin_start
    
    pos = domains[domains['chromosome'] == chromosome][['start', 'end']]
    
    for domain_start, domain_end in pos.values:
        
        #if domain_start < start or domain_end > end: #FIXME?
        
        #    continue
        
        istart = domain_start // binsize
            
        iend = domain_end // binsize

        istart = max(0, istart - bin_start)
        
        iend = min(length, iend - bin_start)

        region = numpy.arange(istart, iend + 1)

        ax.fill_between(region, 
                        region * 0. + ymin, 
                        region * 0. + ymax,
                        color = color,
                        linewidth = 0.,
                        alpha = alpha,
                        clip_on = False)        

    return ax
    
    
def add_scalebar(ax,
                 n_mb,
                 binsize,
                 yfrac = -0.1,
                 style = '-',
                 alpha = 0.9,
                 color = 'black',
                 fontsize = 10,
                 fontweight = 'bold',
                 linewidth = 1.1,
                 bbox_pad = 3):

    #number of bins corresponding to given size
    n_bins = (n_mb * 1e6) / binsize

    #scale bar length (x-axis fraction)
    xfrac = n_bins / ax.get_xlim()[1]
    
    #scale bar position (x-axis fraction)
    pos = 1. - xfrac
    
    #add scale bar
    ax.annotate('', 
                xy = (pos, yfrac), 
                xytext = (1.0, yfrac),                
                xycoords = ax.transAxes,
                textcoords = ax.transAxes,
                annotation_clip = False,
                arrowprops = {'alpha' : alpha, 
                              'color' : color,
                              'arrowstyle' : style,
                              'linewidth' : linewidth})

    #scale bar text position (x-axis fraction)
    tpos = 1. - xfrac / 2.
    
    #add label
    ax.text(tpos, yfrac, 
            '$\mathregular{%s\/Mb}$' % n_mb,
            transform = ax.transAxes,         
            ha = 'center',
            va = 'center',
            alpha = alpha,
            color = color,
            fontsize = fontsize,
            fontweight = fontweight,
            bbox = {'pad' : bbox_pad,
                    'edgecolor' : 'white',
                    'facecolor' : 'white'})

    return ax


def plot_domains(ax, 
                 states,
                 fill_color,
                 alpha = 1.):
    
    #x-axis range
    length = len(states)
    
    #x-axis values (bin positions)
    pos = numpy.arange(length)
    
    #only enriched positions (non-zero)
    enriched = states.replace(0, numpy.nan)
    
    #fill area of enriched regions
    ax.fill_between(pos, 
                    enriched * 0,
                    enriched * 1,
                    color = fill_color,
                    alpha = alpha,
                    linewidth = 0.)

    #set x-axis bounds
    ax.set_xlim(0, length)

    #set y-axis bounds to fill area
    ax.set_ylim(0, 1)

    #hide x-axis ticklabels
    ax.set_xticklabels([])
    
    #hide y-axis ticklabels
    ax.set_yticklabels([])

    for loc in ax.spines:
        #hide spines
        ax.spines[loc].set_visible(False)

    return ax


def distplot(ax,
             values,
             hist_bins,
             hist_color,
             linewidth = 1.2,             
             line_color = greyish,
             line_alpha = 1.,
             dashes = (3.0, 1.5),
             vertical = False,
             zline = -10):

    #plot histogram and kernel density
    sns.distplot(values.dropna(),
                 ax = ax,
                 hist = True,                 
                 bins = hist_bins,
                 color = hist_color,
                 hist_kws = {'edgecolor' : 'white'},
                 vertical = vertical)

    if not vertical:

        #add vertical line
        ax.axvline(linestyle = '--',
                   linewidth = linewidth,
                   color = line_color,
                   alpha = line_alpha,
                   dashes = dashes,
                   zorder = zline)

        #hide x-axis label
        ax.set_xlabel('')

        #hide y-axis ticklabels
        ax.set_yticks([])

        for loc in ['top', 'left', 'right']:
            #hide spines
            ax.spines[loc].set_visible(False)
    
    else:

        #add vertical line
        ax.axhline(linestyle = '--',
                   linewidth = linewidth,
                   color = line_color,
                   alpha = line_alpha,
                   dashes = dashes,
                   zorder = zline)
    
        #hide x-axis label
        ax.set_ylabel('')

        #hide y-axis ticklabels
        ax.set_xticks([])

        for loc in ['top', 'bottom', 'right']:
            #hide spines
            ax.spines[loc].set_visible(False)

    return ax


def kdeplot(ax,
            sample1,
            sample2,
            data,
            cmap,
            n_levels = 10,
            shade = False,
            annotate = True,
            corrfunc = pearsonr,                   
            xytext = (.05, .95),
            fontsize = 10.,
            zorder = 1,
            alpha = 1.,       
            ha = 'left',
            va = 'center',
            bw = 'scott'):

    data = data[[sample1, sample2]].dropna()
    
    if annotate:
    
        #compute correlation coefficient
        r, p = corrfunc(data[sample1], 
                        data[sample2])
        
        #annotate coefficient on plot
        ax.annotate(s = 'r = %.2f' % r, 
                    xy = xytext, 
                    xycoords = ax.transAxes, 
                    ha = ha, 
                    va = va, 
                    fontsize = fontsize)
        
    #plot bivariate kernel density
    sns.kdeplot(data[sample1], 
                data[sample2], 
                ax = ax, 
                bw = bw,
                cmap = cmap,
                shade = shade,
                shade_lowest = False,
                n_levels = n_levels,
                zorder = zorder)
    
    if shade:
        #add contour line if shaded 
        sns.kdeplot(data[sample1], 
                    data[sample2], 
                    ax = ax, 
                    cmap = cmap,
                    shade = False,
                    n_levels = n_levels,
                    zorder = zorder + n_levels,
                    alpha = alpha)
    
    for loc in ['top', 'right']:
        #hide spines
        ax.spines[loc].set_visible(False)

    return ax


