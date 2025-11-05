"""
This script implements my custom matplotlib settings inspired by
https://github.com/fnemina/latexify and adadpted from
mp3_nozzle_experiment_18032022.ipynb.
"""

import os
import matplotlib as mpl
import matplotlib.pyplot  as plt
from cycler import cycler
from collections import OrderedDict
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

fontsize = 17
hwidth = 1.25
length = 5
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

linestyles_dict = {
    'solid': (0, ()),
    'loosely dotted': (0, (1, 10)),
    'dotted': (0, (1, 5)),
    'densely dotted': (0, (1, 1)),
    'loosely dashed': (0, (5, 10)),
    'dashed': (0, (5, 5)),
    'densely dashed': (0, (5, 1)),
    'loosely dashdotted': (0, (3, 10, 1, 10)),
    'dashdotted': (0, (3, 5, 1, 5)),
    'densely dashdotted': (0, (3, 1, 1, 1)),
    'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
    'dashdotdotted': (0, (3, 5, 1, 5, 1, 5)),
    'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1)),
}

params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'savefig.dpi': 600,
    'axes.labelsize': fontsize,
    'axes.titlesize': fontsize,
    'font.size': fontsize,
    'legend.fontsize': fontsize,
    'legend.title_fontsize': fontsize,
    'legend.columnspacing': 0.5, 
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'font.family': 'serif',
    # 'font.serif': ['Computer Modern Roman'],
    'axes.linewidth': hwidth,
    'xtick.major.width': hwidth,
    'ytick.major.width': hwidth,
    'ytick.major.size': length,
    'xtick.major.size': length,
    'xtick.major.pad': 10,  # distance between x-tick and label (in points)
    'ytick.major.pad': 10,  # distance between y-tick and label (in points)
    'xtick.direction': "in",
    'ytick.direction': "in",
    'xtick.top': True,
    'xtick.bottom': True,
    'ytick.left': True,
    'ytick.right': True,
    'grid.alpha': 0.6,
    'grid.linewidth': 0.7 * hwidth,
    'legend.edgecolor': '0',
    'legend.labelspacing': 0.3,
    'axes.edgecolor': 'black',
    'legend.fancybox': False,
    'patch.linewidth': hwidth,
    'legend.facecolor': 'white',
    'legend.framealpha': 1,
    'lines.markersize': 8,
    'lines.markeredgewidth': 1.5,
    'savefig.bbox': 'tight',
    'figure.autolayout': True,
}

mpl.rcParams.update(params)

# Matlab default colours
blue = (0, 0.4470, 0.7410); orange = (0.8500, 0.3250, 0.0980); yellow = (0.9290, 0.6940, 0.1250);
purple = (0.4940, 0.1840, 0.5560); green = (0.4660, 0.6740, 0.1880); lblue = (0.3010, 0.7450, 0.9330);
red = (0.6350, 0.0780, 0.1840); black = (0, 0, 0); dgrey = (0.25, 0.25, 0.25);
grey = (0.5, 0.5, 0.5); lgrey = (0.75, 0.75, 0.75);

# Whittle colours
rob_blue_rgb = '#0072BD' # (0, 114, 189) # shade of blue of Rob's PowerPoint template
whittle_blue_rgb = '#0BACD7' # (11, 172, 215) # shade of blue of Whittle logo

# Colour palette from https://ons-design.notion.site/Colour-335407345de94442b2adccbaa0b0b6e6
plt.rcParams['axes.prop_cycle'] = cycler('color', ['#206095', '#a8bd3a', '#871a5b', '#f66068', '#05341A', '#27a0cc', '#003c57', '#22d0b6', '#746cb1', '#A09FA0'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Add margins around plot
def add_margin(ax, m=0.05):
    for a, s in [(ax.get_xlim, ax.set_xlim), (ax.get_ylim, ax.set_ylim)]:
        lo, hi = a(); r = hi - lo; s(lo - m*r, hi + m*r)

 