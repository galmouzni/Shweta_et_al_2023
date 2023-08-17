#!/usr/bin/env python

import re
import numpy
import pandas
from sklearn.utils import resample
from scipy.stats import pearsonr, spearmanr

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar

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

#functions

def bootstrap_values(matrix,
                     function = numpy.nanmean,
                     confint = 95., 
                     n_boot = 100):
    
    n_cols = matrix.shape[-1]

    sampled_values = numpy.zeros((n_boot, n_cols))
    
    for i in range(n_boot):
        
        resampled = resample(matrix)
        
        sampled_values[i] += function(resampled, 0)
    
    alpha = 100 - confint
    
    lw = numpy.percentile(sampled_values, alpha / 2, 0)
    
    up = numpy.percentile(sampled_values, confint + alpha / 2, 0)
    
    return lw, up
    

def lineplot(ax,
             matrix,
             function = numpy.nanmean,
             confint = 95.,
             n_boot = 100,
             fill_color = blueish,
             line_color = blueish,
             fill_alpha = 0.5,
             line_alpha = 1.0,
             draw_hline = True,
             linewidth = 1.0,
             linestyle = '-',
             zorder = 0):
    
    #x-axis range    
    length = matrix.shape[-1]         

    #set x-axis bounds
    ax.set_xlim(0, length)

    #x-axis range    
    x = numpy.arange(length)
    
    y = function(matrix, 0)
    
    ymin, ymax = bootstrap_values(matrix,
                                  function,
                                  confint, 
                                  n_boot)
        
    ax.plot(x, y, 
            c = line_color, 
            alpha = line_alpha,
            linewidth = linewidth,
            linestyle = linestyle,
            zorder = zorder)

    ax.fill_between(x, 
                    ymin,
                    ymax,
                    linewidth = 0.,
                    color = fill_color,
                    alpha = fill_alpha, 
                    zorder = zorder - 1) 

    if draw_hline:
    
        ax.axhline(linestyle = '--', 
                   linewidth = 0.5,
                   color = 'black',
                   zorder = zorder - 2,
                   alpha = .75, 
                   dashes = (2.5, 2.5))

    #hide x-axis ticks
    ax.set_xticks([])

    #hide y-axis ticks
    ax.set_yticks([])
         
    return ax
    
    
def heatmap(ax,
            matrix,
            vmin = None,
            vmax = None,
            alpha = 1.0,
            cmap = 'RdYlBu_r',
            shading = 'gouraud'):

    if vmin is None:
        #set vmin to 5th percentile if not specified 
        vmin = numpy.nanpercentile(matrix, 5.)
    
    if vmax is None:
        #set vmax to 95th percentile if not specified 
        vmax = numpy.nanpercentile(matrix, 95.)
    
    #plot heatmap
    ax.pcolormesh(matrix,
                  cmap = cmap,
                  alpha = alpha,
                  shading = shading,
                  vmin = vmin,
                  vmax = vmax)
    
    #hide x-axis ticks
    ax.set_xticks([])

    #hide y-axis ticks
    ax.set_yticks([])
    
    return ax
 

def heatmap_im(ax,
               matrix,
               vmin = None,
               vmax = None,
               alpha = 1.0,
               cmap = 'RdYlBu_r',
               interpolation = 'bilinear'):

    if vmin is None:
        #set vmin to 5th percentile if not specified 
        vmin = numpy.nanpercentile(matrix, 5.)
    
    if vmax is None:
        #set vmax to 95th percentile if not specified 
        vmax = numpy.nanpercentile(matrix, 95.)
    
    #plot heatmap
    ax.imshow(matrix[::-1],
              vmin = vmin,
              vmax = vmax,
              cmap = cmap,
              alpha = alpha,
              aspect = 'auto',
              interpolation = interpolation)
    
    #hide x-axis ticks
    ax.set_xticks([])

    #hide y-axis ticks
    ax.set_yticks([])
    
    return ax
 
 
def annotate_line(ax,
                  label1 = 'center', 
                  label2_left = '',
                  label2_right = '',
                  xfrac1 = 0.5,
                  yfrac1 = -0.1,
                  yfrac2 = -0.2,
                  color = 'black',
                  fontsize1 = 10,
                  fontsize2 = 9,
                  fontweight1 = 'medium',
                  fontweight2 = 'medium',
                  arrowprops = {'color' : 'black',
                                'linewidth' : 0.5,
                                'arrowstyle' : '-'},
                  bbox_pad = 1):

    #add scale bar
    ax.annotate('', 
                xy = (0.0, yfrac1), 
                xytext = (1.0, yfrac1),                
                xycoords = ax.transAxes,
                textcoords = ax.transAxes,
                annotation_clip = False,
                arrowprops = arrowprops)
                
    #add label
    ax.text(xfrac1, 
            yfrac1, 
            label1,
            transform = ax.transAxes,         
            ha = 'center',
            va = 'center',
            color = color,
            fontsize = fontsize1,
            fontweight = fontweight1,
            bbox = {'pad' : bbox_pad,
                    'edgecolor' : 'white',
                    'facecolor' : 'white'})
    #add label1
    ax.text(0., yfrac2, 
            label2_left,
            transform = ax.transAxes,         
            ha = 'left',
            va = 'center',
            color = color,
            fontsize = fontsize2,
            fontweight = fontweight2)
    
    #add label2
    ax.text(1., yfrac2, 
            label2_right,
            transform = ax.transAxes,         
            ha = 'right',
            va = 'center',
            color = color,
            fontsize = fontsize2,
            fontweight = fontweight2)

    return ax


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
  
    norm = colors.Normalize(vmin, vmax)
     
    cbar = colorbar.ColorbarBase(ax, 
                                 cmap,
                                 norm,
                                 orientation = orientation)
    
    ax.tick_params(pad = label_pad,
                   length = tick_length)
       
    cbar.set_ticks(ticks)
    
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
    return ax, cbar
  

