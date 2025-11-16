import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from cycler import cycler
from matplotlib import gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
from shapely import Polygon

from matplotlib_custom_settings import *

# xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'point_load_study_results_regional_spanfrac.xlsx'))
# xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'point_load_study_results_regional_spanfrac_fuse.xlsx'))
xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'point_load_study_results_regional_spanfrac_wing_new.xlsx'))
# xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'mass', 'point_load_study_results_regional_spanfrac.xlsx'))

sheet_dict = {}
for sheet_name in xlsx.sheet_names:
    sheet_dict[sheet_name] = pd.read_excel(xlsx, sheet_name=sheet_name)

#%%

df = sheet_dict['results']

fcs_loc_unique = np.sort(df["fcs_loc"].unique())
span_loc_unique = np.sort(df["span_loc"].unique())
wing_frac_unique = np.sort(df["wing_frac"].unique())
sigma_fc_unique = np.sort(df["sigma_fc"].unique())

# Create 2D coordinate grids
A, B, C, D = np.meshgrid(fcs_loc_unique, span_loc_unique, wing_frac_unique, sigma_fc_unique, indexing="ij")

# Create shape tuple for reshaping
shape = (
    len(fcs_loc_unique),
    len(span_loc_unique),
    len(wing_frac_unique),
    len(sigma_fc_unique),
)

# Create index mapping for each dimension
fcs_loc_idx = {v: i for i, v in enumerate(fcs_loc_unique)}
span_loc_idx = {v: i for i, v in enumerate(span_loc_unique)}
wing_frac_idx = {v: i for i, v in enumerate(wing_frac_unique)}
sigma_fc_idx = {v: i for i, v in enumerate(sigma_fc_unique)}

# Initialize empty arrays
CDS_grid = np.full(shape, np.nan)
CDS_ref_grid = np.full(shape, np.nan)
WMTO_grid = np.full(shape, np.nan)
WMTO_ref_grid = np.full(shape, np.nan)
Wpointload_fuselage_grid = np.full(shape, np.nan)
fcs_loc_grid = np.full(shape, np.nan)
span_loc_grid = np.full(shape, np.nan)
wing_frac_grid = np.full(shape, np.nan)
sigma_fc_grid = np.full(shape, np.nan)

# Populate arrays
for _, row in df.iterrows():
    i = fcs_loc_idx[row["fcs_loc"]]
    j = span_loc_idx[row["span_loc"]]
    k = wing_frac_idx[row["wing_frac"]]
    l = sigma_fc_idx[row["sigma_fc"]]

    CDS_grid[i, j, k, l] = row["CDS"]
    CDS_ref_grid[i, j, k, l] = row["CDS_ref"]
    WMTO_grid[i, j, k, l] = row["WMTO"]
    WMTO_ref_grid[i, j, k, l] = row["WMTO_ref"]
    Wpointload_fuselage_grid[i, j, k, l] = row["Wpointload_fuselage"]
    
    fcs_loc_grid[i, j, k, l] = row["fcs_loc"]
    span_loc_grid[i, j, k, l] = row["span_loc"]
    wing_frac_grid[i, j, k, l] = row["wing_frac"]
    sigma_fc_grid[i, j, k, l] = row["sigma_fc"]
    
    
#%%

# idx_fcs_loc_des = -1  # inactive
# idx_span_loc_des = 1  # inactive
# idx_wing_frac_des = -1
# idx_sigma_fc_des = 0
    
# fig, ax = plt.subplots()
    
# ax.scatter(
#     fcs_loc_grid[:, :, 2, -1], span_loc_grid[:, :, 2, -1],
#     c=WMTO_grid[:, :, 2, -1],
#     marker='o',
#     # color='black',
#     zorder=20,
# )

# plt.show()
    
    
# sys.exit()

#%%

idx_fcs_loc_des = -1  # inactive
# idx_span_loc_des = 1  # inactive
idx_span_loc_des = 0  # inactive
idx_wing_frac_des = -1
idx_sigma_fc_des = 0

fig, ax = plt.subplots(figsize=(8, 6))

"""

ax.scatter(
    WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    marker='o',
    color='black',
    zorder=20,
)

ax.plot(
    WMTO_grid[:, idx_span_loc_des, 0, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[:, idx_span_loc_des, 0, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[0], marker='.', zorder=100,
)

ax.plot(
    WMTO_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[1], marker='.',
)

ax.plot(
    WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[2], marker='.',
)

ax.plot(
    WMTO_grid[0, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[0, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[2], marker='.',
)
ax.plot(
    WMTO_grid[1, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[1, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[2], marker='.',
)
ax.plot(
    WMTO_grid[2, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[2, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[2], marker='.',
)
ax.plot(
    WMTO_grid[3, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[3, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[2], marker='.',
)

ax.plot(
    WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[3], marker='.',
)

# cvhull_points = np.array([
#     WMTO_grid[:, :-1, :, idx_sigma_fc_des].flatten() / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des].flatten(),
#     CDS_grid[:, :-1, :, idx_sigma_fc_des].flatten() / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax.fill(
#     x_outline_polygon, y_outline_polygon, color='grey', alpha=0.1, zorder=-1,
# )
"""
#####

idx_sigma_fc_des = -1

# ax.scatter(
#     WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     marker='o',
#     color='black',
#     zorder=20,
# )

# ax.plot(
#     WMTO_grid[:, idx_span_loc_des, 0, idx_sigma_fc_des],# / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[:, idx_span_loc_des, 0, idx_sigma_fc_des],# / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[0], marker='.', zorder=100,
# )

ax.plot(
    WMTO_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des],# / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des],# / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[1], marker='.',
)

# ax.plot(
#     WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[2], marker='.',
# )

# ax.plot(
#     WMTO_grid[0, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[0, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[2], marker='.',
# )
# ax.plot(
#     WMTO_grid[1, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[1, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[2], marker='.',
# )
# ax.plot(
#     WMTO_grid[2, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[2, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[2], marker='.',
# )
# ax.plot(
#     WMTO_grid[3, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[3, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[2], marker='.',
# )

# ax.plot(
#     WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     CDS_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
#     color=colors[3], marker='.',
# )

# cvhull_points = np.array([
#     WMTO_grid[:, :-1, :, -1].flatten() / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des].flatten(),
#     CDS_grid[:, :-1, :, -1].flatten() / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax.fill(
#     x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
# )

#####
"""
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Create zoomed inset
axins = inset_axes(ax, width="100%", height="100%", loc='lower left',
                   bbox_to_anchor=(1.1, 0.725, 0.4, 0.3),
                   bbox_transform=ax.transData, borderpad=0)

# Plot same data in inset
for line in ax.get_lines():
    axins.plot(
        line.get_xdata(), line.get_ydata(),
        color=line.get_color(),
        linestyle=line.get_linestyle(),
        marker=line.get_marker(),
        alpha=line.get_alpha() if line.get_alpha() is not None else 1.0,
        zorder=line.get_zorder(),
    )
    
# from copy import copy
# for patch in ax.patches:
#     p = copy(patch)
#     axins.add_patch(p)

from matplotlib.patches import Polygon

# Replot all patch artists (alpha-shape polygons) safely
for patch in ax.patches:
    if isinstance(patch, Polygon):
        new_patch = Polygon(
            patch.get_xy(),
            closed=True,
            facecolor=patch.get_facecolor(),
            edgecolor=patch.get_edgecolor(),
            alpha=patch.get_alpha(),
            zorder=patch.get_zorder(),
        )
        axins.add_patch(new_patch)

for coll in ax.collections:
    # offsets are Nx2 data coordinates
    offsets = coll.get_offsets()
    xs = offsets[:, 0]
    ys = offsets[:, 1]
    # sizes, facecolors, edgecolors, alpha
    sizes = coll.get_sizes()
    # If sizes is empty or single value, scatter will handle it correctly
    facecolors = coll.get_facecolors()
    edgecolors = coll.get_edgecolors()
    alpha = coll.get_alpha()
    # try to preserve marker shape if available
    marker_path = coll.get_paths()[0]  # Path for marker

    # call scatter on inset: preserve face/edgecolors and sizes when possible
    axins.scatter(xs, ys, s=sizes, facecolors=facecolors, edgecolors=edgecolors,
                  marker=marker_path, alpha=alpha, zorder=coll.get_zorder())

# Zoomed limits
axins.set_xlim(1.22, 1.43)
axins.set_ylim(1.11, 1.33)
axins.set_xticks([])
axins.set_yticks([])

# Connect inset with main plot
patch, connector1, connector2 = mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="black", alpha=0.2)

from matplotlib.patches import Rectangle
bbox = axins.get_position()  # get inset box in figure coordinates
bg = Rectangle(
    (0, 0), 1, 1,
    transform=axins.transAxes,  # full inset area
    facecolor="white",
    edgecolor="none",
    zorder=-5  # put *behind all data* but *above connectors*
)
axins.add_patch(bg)
bg.set_zorder(-5)
patch.set_zorder(5)

# Now adjust the zorder of the connector lines:
for c in [connector1, connector2]:
    c.set_zorder(-100)   # draw behind the boxes
    c.set_alpha(0.2)   # draw behind the boxes
"""
#####

# ax.plot([0.7, 1.6], [0.7, 1.6], color='black', alpha=0.2)
# ax.set_xlim(0.75, 1.5)
# ax.set_ylim(0.7, 1.4)

#####

legend_elements = [
    Line2D([0], [0], color='black', marker='o', label=r'Baseline'),#' $\left( \hat{y}=0.37, \%_\mathrm{wing}=1, \sigma_\mathrm{FC}=2~\mathrm{kW/kg} \right)$'),
    Line2D([0], [0], color=colors[0], label=r'$\hat{x} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[1], label=r'$\hat{y} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[2], label=r'$\%_\mathrm{wing} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[3], label=r'$\sigma_\mathrm{FC} \rightarrow \left[ 2,4 \right]$ kW/kg'),
    # Line2D([0], [0], color='grey', linestyle='solid', label=r'$\hat{x}=c$'),
    # Line2D([0], [0], color='grey', linestyle='dashed', label=r'$N_\mathrm{eng}=c$'),
    # Line2D([0], [0], color='grey', linestyle='dashdot', label=r'$\%_\mathrm{wing}=c$'),
]
leg_main = ax.legend(handles=legend_elements, loc='upper left', frameon=False)
ax.add_artist(leg_main)

ax.set_xlabel("Relative weight (-)", labelpad=10)
ax.set_ylabel("Relative drag (-)", labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
    
# plt.savefig("weight_vs_drag_spanfrac.svg", format='svg')

plt.show()

sys.exit()




















#%%


























































#%%

def add_margin(ax, m=0.05):
    for a, s in [(ax.get_xlim, ax.set_xlim), (ax.get_ylim, ax.set_ylim)]:
        lo, hi = a(); r = hi - lo; s(lo - m*r, hi + m*r)

plt.rcParams['axes.prop_cycle'] = cycler('color', ['#206095', '#a8bd3a', '#871a5b', '#f66068', '#05341A', '#27a0cc', '#003c57', '#22d0b6', '#746cb1', '#A09FA0'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

legend_elements = [
    Line2D([0], [0], color='black', marker='o', label=r'Baseline $\left( \hat{y}=0.37, \%_\mathrm{wing}=1, \sigma_\mathrm{FC}=2~\mathrm{kW/kg} \right)$'),
    Line2D([0], [0], color=colors[0], label=r'$\hat{x} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[1], label=r'$\hat{y} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[2], label=r'$\%_\mathrm{wing} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[3], label=r'$\sigma_\mathrm{FC} \rightarrow \left[ 2,4 \right]$ kW/kg'),
    # Line2D([0], [0], color='grey', linestyle='solid', label=r'$\hat{x}=c$'),
    # Line2D([0], [0], color='grey', linestyle='dashed', label=r'$N_\mathrm{eng}=c$'),
    # Line2D([0], [0], color='grey', linestyle='dashdot', label=r'$\%_\mathrm{wing}=c$'),
]

idx_type = 'start'
# idx_type = 'end'
if idx_type == 'start':
    idx = 0
elif idx_type == 'end':
    idx = -1
    
sigma_fc_idx = 0

fig, ax = plt.subplots(figsize=(8, 6))

# Explicit axis limits needed when not using ax.scatter() due to nan's in some arrays
X = WMTO_grid / WMTO_ref_grid
Y = CDS_grid / CDS_ref_grid
"""
ax.set_xlim(np.nanmin(X), np.nanmax(X))
ax.set_ylim(np.nanmin(Y), np.nanmax(Y))
"""
ax.set_xlim(0.96, 1.3)
ax.set_ylim(0.96, 1.2)

xlim = ax.get_xlim()
ylim = ax.get_ylim()

##
"""
cvhull_points = np.array([
    WMTO_grid.flatten() / WMTO_ref_grid.flatten(),
    CDS_grid.flatten() / CDS_ref_grid.flatten(),
]).T
cvhull = ConvexHull(cvhull_points)
cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
polygon = Polygon(cvhull_vertices)
x_outline_polygon, y_outline_polygon = polygon.exterior.xy
ax.fill(
    x_outline_polygon, y_outline_polygon, color='grey', alpha=0.1, zorder=-1,
)

##

cvhull_points = np.array([
    WMTO_grid[:,:,:,0].flatten() / WMTO_ref_grid[:,:,:,0].flatten(),
    CDS_grid[:,:,:,0].flatten() / CDS_ref_grid[:,:,:,0].flatten(),
]).T
cvhull = ConvexHull(cvhull_points)
cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
polygon = Polygon(cvhull_vertices)
x_outline_polygon, y_outline_polygon = polygon.exterior.xy
ax.fill(
    x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
)
"""
##

cvhull_points = np.array([
    WMTO_grid[:,:,:,-1].flatten() / WMTO_ref_grid[:,:,:,-1].flatten(),
    CDS_grid[:,:,:,-1].flatten() / CDS_ref_grid[:,:,:,-1].flatten(),
]).T
cvhull = ConvexHull(cvhull_points)
cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
polygon = Polygon(cvhull_vertices)
x_outline_polygon, y_outline_polygon = polygon.exterior.xy
ax.fill(
    x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
)

##

# ax.scatter(WMTO_grid / WMTO_ref_grid, CDS_grid / CDS_ref_grid, marker='.')

ax.set_xlabel("Relative weight (-)", labelpad=10)
ax.set_ylabel("Relative drag (-)", labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

"""
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)

ax.set_xlim(left=1)
ax.set_ylim(bottom=1)
"""
## Legend

leg_main = ax.legend(handles=legend_elements, loc='upper left', frameon=False)
ax.add_artist(leg_main)

# Create a proxy artist for the legend
space_patch_1 = mpatches.Patch(color='grey', alpha=0.1, label='Possibility space')
space_patch_2 = mpatches.Patch(color='black', alpha=0.2, label='Low-tech FC design space')
space_patch_3 = mpatches.Patch(color='black', alpha=0.2, label='High-tech FC design space')

space_handles = [space_patch_1, space_patch_2, space_patch_3]
leg_space = ax.legend(handles=space_handles, loc='lower right', frameon=False)

for leg in (leg_main, leg_space):
    leg.set_zorder(100)        # draw legends on top of plot
    leg.get_frame().set_alpha(0)  # transparent background (if frameon=True)

# ax.set_xscale('log')
# ax.set_yscale('log')

#####

idx_fcs_loc_des = -1  # inactive
idx_span_loc_des = 1  # inactive
idx_wing_frac_des = -1
idx_sigma_fc_des = 0

ax.scatter(
    WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, idx_sigma_fc_des],
    marker='o',
    color='black',
    zorder=20,
)

for idx_fcs_loc in range(len(fcs_loc_unique)):
    # ax.plot(
    #     WMTO_grid[idx_fcs_loc, idx_span_loc_des, [0, -1], idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc, idx_span_loc_des, [0, -1], idx_sigma_fc_des],
    #     CDS_grid[idx_fcs_loc, idx_span_loc_des, [0, -1], idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, [0, -1], idx_sigma_fc_des],
    #     color=colors[0],
    #     marker='.',
    # )
    ax.plot(
        WMTO_grid[idx_fcs_loc, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc, idx_span_loc_des, :, idx_sigma_fc_des],
        CDS_grid[idx_fcs_loc, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des],
        color=colors[0],
    )
    
    ax.scatter(
        WMTO_grid[idx_fcs_loc, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc, idx_span_loc_des, :, idx_sigma_fc_des],
        CDS_grid[idx_fcs_loc, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des],
        color=colors[2],
        marker='.',
        zorder=10,
    )
    
# for idx_span_loc in range(len(span_loc_unique)):
ax.plot(
    WMTO_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des],
    CDS_grid[idx_fcs_loc_des, :, idx_wing_frac_des, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des_des, :, idx_wing_frac_des, idx_sigma_fc_des],
    color=colors[1],
    marker='.',
)
    
# # for idx_wing_frac in range(len(wing_frac_unique)):
# ax.scatter(
#     WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des],
#     CDS_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, :, idx_sigma_fc_des],
#     color=colors[2],
#     marker='.',
# )
"""
# for idx_sigma_fc in range(len(sigma_fc_unique)):
ax.plot(
    WMTO_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :] / WMTO_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :],
    CDS_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :] / CDS_ref_grid[idx_fcs_loc_des, idx_span_loc_des, idx_wing_frac_des, :],
    color=colors[3],
    marker='.',
)

cvhull_points = np.array([
    WMTO_grid[:,:,:,0].flatten() / WMTO_ref_grid[:,:,:,0].flatten(),
    CDS_grid[:,:,:,0].flatten() / CDS_ref_grid[:,:,:,0].flatten(),
]).T
cvhull = ConvexHull(cvhull_points)
cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
polygon = Polygon(cvhull_vertices)
x_outline_polygon, y_outline_polygon = polygon.exterior.xy
ax.fill(
    x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
)
"""
add_margin(ax, 0.05)
xlim = ax.get_xlim()
ax.set_ylim(bottom=xlim[0])
ax.plot([0, 2], [0, 2], color='black', linestyle='dashed', alpha=0.1)
ax.set_aspect('equal')

# plt.tight_layout()
# plt.savefig('figure.svg', format='svg')
plt.show()

sys.exit()

#%%
































#%%

# Z = WMTO_grid[:,:,-1,0] / WMTO_ref_grid[:,:,-1,0]
# Z = WMTO_grid[:,:,0,0] / WMTO_ref_grid[:,:,0,0]
Z = WMTO_grid[:,:,1,0] / WMTO_ref_grid[:,:,1,0]

fig, ax = plt.subplots()

# cf = ax.contourf(A[:,:,-1,0], B[:,:,-1,0], Z)
# cf = ax.contourf(A[:,:,0,0], B[:,:,0,0], Z)
cf = ax.contourf(A[:,:,1,0], B[:,:,1,0], Z)
cbar = fig.colorbar(cf, ax=ax, label="Relative weight (-)")

plt.show()


#%%

# fig, ax = plt.subplots()

# for i in range(4):
#     for j in range(4):
#         ax.scatter(D[:,:,1,j], WMTO_grid[:,:,1,j] / WMTO_ref_grid[:,:,1,j], color=colors[j])

# plt.show()


#%%

fig, ax = plt.subplots()

X = D[:,0,0,:]
Y = WMTO_grid[:,0,0,:] / WMTO_ref_grid[:,0,0,:]

ax.scatter(X, Y, color=colors[0])

X = D[0,:,0,:]
Y = WMTO_grid[0,:,0,:] / WMTO_ref_grid[0,:,0,:]

ax.scatter(X, Y, color=colors[1])

X = D[0,0,:,:]
Y = WMTO_grid[0,0,:,:] / WMTO_ref_grid[0,0,:,:]

ax.scatter(X, Y, color=colors[2])

plt.show()



