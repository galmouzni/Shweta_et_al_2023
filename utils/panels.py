#!/usr/bin/env python

from matplotlib import pyplot

class Panel():

    def __init__(self,
                 figure,
                 panel_left,
                 panel_bottom,
                 panel_width,
                 panel_height):

        #panel left margin (inches)
        self.left = panel_left

        #panel bottom margin (inches)
        self.bottom = panel_bottom

        #panel width (inches)
        self.width = panel_width
        
        #panel height (inches)
        self.height = panel_height

        #figure width (inches)
        self.fig_width = figure.width

        #figure height (inches)
        self.fig_height = figure.height

        #determine position (as fraction of figure inches)
        position = self.position()
        
        #add panel as new axis at given position
        self.axis = figure.plot.add_axes(position)

    @property
    def left_position(self):
        #left position relative to figure width
        return self.left / self.fig_width

    @property
    def bottom_position(self):    
        #bottom position relative to figure height
        return self.bottom / self.fig_height
  
    @property
    def relative_width(self):    
        #width relative to figure height
        return self.width / self.fig_width  
    
    @property
    def relative_height(self): 
        #height relative to figure height
        return self.height / self.fig_height
    
    def position(self):
        #panel position
        return (self.left_position,
                self.bottom_position,
                self.relative_width,
                self.relative_height)
    
    def update_position(self,
                        new_width,
                        new_height):

        #update figure width
        self.fig_width = new_width
        
        #update figure height
        self.fig_height = new_height
        
        #determine position
        position = self.position()
        
        #update panel axis position
        self.axis.set_position(self.position())


class MultiPanel():

    def __init__(self,
                 top = 0.5,
                 bottom = 0.5,
                 left = 0.5,
                 right = 0.5,
                 hspace = 0.1,
                 vspace = 0.1,
                 panel_width = 1.0,
                 panel_height = 1.0):
        
        #figure margins (inches)
        self.top = top
        self.bottom = bottom
        self.left = left
        self.right = right

        #figure width (inches)
        self.width = left + panel_width + right
        
        #figure height (inches)
        self.height = bottom + panel_height + top

        #horizontal space between panels (inches)
        self.hspace = hspace

        #vertical space between panels (inches)
        self.vspace = vspace

        #default panel width (inches)
        self.panel_width = panel_width

        #default panel height (inches)
        self.panel_height = panel_height

        #create figure
        self.plot = pyplot.figure()

        #initialize first panel
        panel = Panel(self,
                      left,
                      bottom,
                      panel_width,
                      panel_height) 

        #list of panels in figure
        self.panels = [panel]

        #list of axes in figure
        self.axes = [panel.axis]
        
        #update figure size
        self.update_figure()

    def update_figure(self):

        #update figure size
        self.plot.set_size_inches(self.width, 
                                  self.height)

    def add_panel(self,
                  panel_left,
                  panel_bottom,
                  panel_width,
                  panel_height,
                  below = True):

        #create new panel
        new_panel = Panel(self,
                          panel_left,
                          panel_bottom,
                          panel_width,
                          panel_height)        

        #add new panel to list
        self.panels.append(new_panel)
        
        #add new axis to list
        self.axes.append(new_panel.axis)
        
        #if new panel below bottom margin
        if below and panel_bottom < self.bottom:

            for panel in self.panels:
                #shift panels upwards
                panel.bottom += panel_height + self.hspace
            
            #increase figure height
            self.height += panel_height + self.hspace
        
        #if new panel larger than figure
        if panel_left + panel_width + self.right > self.width:
            #increase figure width
            self.width += panel_width + self.vspace
        
        #update figure size
        self.update_figure()

        #update panel position relative to new figure size
        for panel in self.panels:

            panel.update_position(self.width,
                                  self.height)
        
    def add_below(self,
                  panel,
                  space_frac = 1.0,
                  width_frac = 1.0,
                  height_frac = 1.0):

        panel_width = self.panel_width * width_frac

        panel_height = self.panel_height * height_frac

        #left margin same as given panel
        panel_left = panel.left
        
        #bottom margin (horizontal space plus panel height below given panel)
        panel_bottom = panel.bottom - panel_height - self.hspace * space_frac
        
        self.add_panel(panel_left,
                       panel_bottom,
                       panel_width,
                       panel_height)
        
        return self.panels[-1]
        
    def add_after(self,
                  panel,
                  space_frac = 1.0,                  
                  width_frac = 1.0,
                  height_frac = 1.0):

        panel_width = self.panel_width * width_frac

        panel_height = self.panel_height * height_frac

        #left margin (vertical space after given panel)
        panel_left = panel.left + panel.width + self.vspace * space_frac #FIXME this is problematic if panel is added to the right of a panel added below
        
        #bottom margin same as given panel
        panel_bottom = panel.bottom
 
        self.add_panel(panel_left,
                       panel_bottom,
                       panel_width,
                       panel_height,
                       below = False)
                                              
        return self.panels[-1]

