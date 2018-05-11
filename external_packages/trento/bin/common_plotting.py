import matplotlib
import matplotlib.pyplot as plt

import ConfigParser

import argparse

from itertools import cycle

default_colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

default_line_sytles = [
        #dash: 20 points on, 20 points off
        [20, 10],
        #small-dash: 10 points on, 5 points off
        [10, 5],
        #dash-dotted line: 8 on, 4 off, 2 on, 4 off
        [8, 4, 2, 4],
        #dash-dot-dot line
        [8, 4, 2, 4, 2, 4],
        #dot-dash-dash line
        [2, 4, 8, 4, 8, 4],
        #dotted line: 2 points on, 2 points off
        [2, 2]
        ]


smash_style_default = {'backend' : 'pdf',
        'axes.labelsize': 40,
        'axes.titlesize': 40,
        'axes.linewidth' : 2.0,
        'axes.unicode_minus': True,
        'axes.color_cycle' : default_colors,

        'legend.frameon': False,
        'legend.fontsize': 30,
        'legend.numpoints': 1,

        'lines.markeredgewidth'  : 5.0,
        'lines.linewidth' : 5.0,
        'lines.markersize'  : 10,
        'lines.markeredgewidth'  : 2.5,

        'mathtext.fontset': 'stix',
        #'text.usetex': True,  - crashes
        #'font.family':'serif',
        'font.serif':['Computer Modern Roman', 'CMU Serif'],
        'font.monospace': ['Computer Modern Typewriter', 'CMU Typewriter Text'],
        'font.size' : 30,
        'font.family':'sans-serif',
        'mathtext.default':'rm',

        'xtick.labelsize': 30,
        'xtick.major.size':10,
        'xtick.minor.size':5,
        'xtick.major.pad':15,
        'xtick.major.width' : 2.0,

        'ytick.labelsize': 30,
        'ytick.major.size':10,
        'ytick.minor.size':5,
        'ytick.major.pad':15,
        'ytick.major.width' : 2.0,

        'figure.figsize': (15., 8.)}


class SmashStyle(argparse.Namespace):
    '''Smash style plot template for matplotlib

    The smash style provides proper figsize, fontsize, fonttype, linewidth,
    color cycles, ticksize and padding and legend style.

    Example usage:
    If one want to use the basic styles except linesytles, increased title
    and plot space, minorticks, one can import the template file in the
    beginning of python script:

          from common_plotting import smash_style

    or simply:

          import common_plotting

    Otherwise, call the set function before plt.legend() if there is one
    or before plt.show() or plt.savefig() if there is no plt.legend() like:

          smash_style.set()
    '''


    def __init__(self):
        self.params = smash_style_default
        plt.rcParams.update(self.params)

    #get all figs in the current plot
    def get_all_figs(self):
        fignums = plt.get_fignums()
        figs = [plt.figure(i) for i in fignums]

        return figs

    #get all subplots in one fig
    def get_all_subplots(self, fig):
        return fig.get_axes()

    #get all lines in one subplots
    def get_all_lines(self, axes):
        children = axes.get_children()
        lines = [l for l in children if type(l)==matplotlib.lines.Line2D]
        return lines

    #get all legends in one ax
    def update_legends(self, ax):
        children = ax.get_children()
        legends = [l for l in children if type(l)==matplotlib.legend.Legend]

        for legend in legends:
            title = legend.get_title().get_text()
            handles, labels = ax.get_legend_handles_labels()
            if title == 'None':
                ax.legend(handles, labels)
            else:
                ax.legend(handles, labels, title=title)

    #apply the smash style to the current script
    def set(self, line_styles=True, title_padding=1.03, minorticks_on=True):
        '''Args:
             line_styles (bool, optional): True to use smash line style cycles.
                 Defaults to True.
             title_padding (float, optional): To increase padding between title
                 and plot. Defaults to 1.03, bigger title_padding to move title
                 upper.
             minorticks_on: True to turn on minor ticks. Defaults to True.
                 False to switch off minorticks'''
        figs = self.get_all_figs()
        for fig in figs:
            axes = self.get_all_subplots(fig)
            for ax in axes:
                if minorticks_on:
                    ax.minorticks_on()

                if title_padding:
                    #increase the title padding
                    title_text = ax.get_title()
                    ax.set_title(title_text, y=title_padding)

                if line_styles:
                    lines = self.get_all_lines(ax)

                    line_style_cycle = cycle(default_line_sytles)
                    for i, line in enumerate(lines):
                        #skip if line.markers=None
                        if line.get_linestyle() == 'None':
                            continue
                        #for i==0, use solid line
                        if i > 0:
                            plt.setp(line, dashes=next(line_style_cycle))

                    self.update_legends(ax)



smash_style = SmashStyle()
