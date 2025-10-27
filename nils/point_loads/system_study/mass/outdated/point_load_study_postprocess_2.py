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

xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'exported_study_results_corr.xlsx'))

sheet_dict = {}
for sheet_name in xlsx.sheet_names:
    sheet_dict[sheet_name] = pd.read_excel(xlsx, sheet_name=sheet_name)

#%%

df = sheet_dict['results']

fcs_loc_unique = np.sort(df["fcs_loc"].unique())
N_eng_unique = np.sort(df["N_eng"].unique())
wing_frac_unique = np.sort(df["wing_frac"].unique())
sigma_fc_unique = np.sort(df["sigma_fc"].unique())

# Create 2D coordinate grids
A, B, C, D = np.meshgrid(fcs_loc_unique, N_eng_unique, wing_frac_unique, sigma_fc_unique, indexing="ij")

# Create shape tuple for reshaping
shape = (
    len(fcs_loc_unique),
    len(N_eng_unique),
    len(wing_frac_unique),
    len(sigma_fc_unique),
)

# Create index mapping for each dimension
fcs_loc_idx = {v: i for i, v in enumerate(fcs_loc_unique)}
N_eng_idx = {v: i for i, v in enumerate(N_eng_unique)}
wing_frac_idx = {v: i for i, v in enumerate(wing_frac_unique)}
sigma_fc_idx = {v: i for i, v in enumerate(sigma_fc_unique)}

# Initialize empty arrays
CDS_grid = np.full(shape, np.nan)
CDS_ref_grid = np.full(shape, np.nan)
WMTO_grid = np.full(shape, np.nan)
WMTO_ref_grid = np.full(shape, np.nan)
Wpointload_fuselage_grid = np.full(shape, np.nan)

# Populate arrays
for _, row in df.iterrows():
    i = fcs_loc_idx[row["fcs_loc"]]
    j = N_eng_idx[row["N_eng"]]
    k = wing_frac_idx[row["wing_frac"]]
    l = sigma_fc_idx[row["sigma_fc"]]

    CDS_grid[i, j, k, l] = row["CDS"]
    CDS_ref_grid[i, j, k, l] = row["CDS_ref"]
    WMTO_grid[i, j, k, l] = row["WMTO"]
    WMTO_ref_grid[i, j, k, l] = row["WMTO_ref"]
    Wpointload_fuselage_grid[i, j, k, l] = row["Wpointload_fuselage"]

#%%

def add_margin(ax, m=0.05):
    for a, s in [(ax.get_xlim, ax.set_xlim), (ax.get_ylim, ax.set_ylim)]:
        lo, hi = a(); r = hi - lo; s(lo - m*r, hi + m*r)

plt.rcParams['axes.prop_cycle'] = cycler('color', ['#206095', '#a8bd3a', '#871a5b', '#f66068', '#05341A', '#27a0cc', '#003c57', '#22d0b6', '#746cb1', '#A09FA0'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

legend_elements = [
    Line2D([0], [0], color=colors[0], label=r'$\left( \hat{x}, N_\mathrm{eng} \right)$'),
    Line2D([0], [0], color=colors[1], label=r'$\left( \hat{x}, \%_\mathrm{wing} \right)$'),
    Line2D([0], [0], color=colors[2], label=r'$\left(N_\mathrm{eng}, \%_\mathrm{wing} \right)$'),
    Line2D([0], [0], color='grey', linestyle='solid', label=r'$\hat{x}=c$'),
    Line2D([0], [0], color='grey', linestyle='dashed', label=r'$N_\mathrm{eng}=c$'),
    Line2D([0], [0], color='grey', linestyle='dashdot', label=r'$\%_\mathrm{wing}=c$'),
]

idx_type = 'start'
# idx_type = 'end'
if idx_type == 'start':
    idx = 0
elif idx_type == 'end':
    idx = -1
    
sigma_fc_idx = 0

# fig, ax = plt.subplots(figsize=(8, 6))

fig = plt.figure(figsize=(12,9), constrained_layout=True)
gs = gridspec.GridSpec(
    2, 2, width_ratios=[5,2], height_ratios=[2,5], hspace=0.1, wspace=0.1,
)

ax_main = fig.add_subplot(gs[1,0])
ax_blue = fig.add_subplot(gs[0,0])
ax_green = fig.add_subplot(gs[1,1])
ax_purple = fig.add_subplot(gs[0,1])

for sigma_fc_idx in [0, -1]:

    ## Vary fcs_loc_unique and N_eng_unique
    
    idx = 0
    
    X = WMTO_grid[:,:,idx,sigma_fc_idx] / WMTO_ref_grid[:,:,idx,sigma_fc_idx]
    Y = CDS_grid[:,:,idx,sigma_fc_idx] / CDS_ref_grid[:,:,idx,sigma_fc_idx]
    
    fcs_loc_unique = fcs_loc_unique.astype(float)
    fcs_loc_unique[-1] -= 1e-6
    cs = ax_main.contour(
        X,
        Y,
        A[:,:,idx,sigma_fc_idx],
        levels=fcs_loc_unique,
        colors=colors[0],
        linestyles='solid',
        linewidths=0.8,
        zorder=10,
        extend='max',
    )
    N_eng_unique = N_eng_unique.astype(float)
    N_eng_unique[-1] -= 1e-6
    cs2 = ax_main.contour(
        X,
        Y,
        B[:,:,idx,sigma_fc_idx],
        levels=N_eng_unique,
        colors=colors[0],
        linestyles='dashed',
        linewidths=0.8,
        zorder=10,
        extend='max',
    )
    ax_main.contourf(
        X,
        Y,
        A[:,:,idx,sigma_fc_idx],
        levels=fcs_loc_unique,
        colors=colors[0],
        alpha=0.1,
        zorder=5,  # below your contour lines
    )
    
    # ax_main.scatter(WMTO_grid[:,:,0,0] / WMTO_ref_grid[:,:,0,0], CDS_grid[:,:,0,0] / CDS_ref_grid[:,:,0,0], marker='.')
    
    # Explicit axis limits needed when not using ax_main.scatter() due to nan's in some arrays
    # ax_main.set_xlim(np.nanmin(X), np.nanmax(X))
    # ax_main.set_ylim(np.nanmin(Y), np.nanmax(Y))
    
    ## Vary fcs_loc_unique and wing_frac_idx
    
    idx = 0
    
    X = WMTO_grid[:,idx,:,sigma_fc_idx] / WMTO_ref_grid[:,idx,:,sigma_fc_idx]
    Y = CDS_grid[:,idx,:,sigma_fc_idx] / CDS_ref_grid[:,idx,:,sigma_fc_idx]
    
    fcs_loc_unique = fcs_loc_unique.astype(float)
    fcs_loc_unique[-1] -= 1e-6
    cs = ax_main.contour(
        X,
        Y,
        A[:,idx,:,sigma_fc_idx],
        levels=fcs_loc_unique,
        colors=colors[1],
        linestyles='solid',
        linewidths=0.8,
        zorder=10,
        extend='max',
    )
    wing_frac_unique = wing_frac_unique.astype(float)
    wing_frac_unique[-1] -= 1e-6
    cs2 = ax_main.contour(
        X,
        Y,
        C[:,idx,:,sigma_fc_idx],
        levels=wing_frac_unique,
        colors=colors[1],
        linestyles='dashdot',
        linewidths=0.8,
        zorder=10,
        extend='max',
    )
    # Fill the area covered by the first contour set with red, alpha=0.5
    ax_main.contourf(
        X,
        Y,
        A[:,idx,:,sigma_fc_idx],
        levels=fcs_loc_unique,
        colors=colors[1],
        alpha=0.1,
        zorder=5,  # below your contour lines
    )
    
    ## Vary N_eng_idx and wing_frac_idx
    
    idx = -1
    
    X = WMTO_grid[idx,:,:,sigma_fc_idx] / WMTO_ref_grid[idx,:,:,sigma_fc_idx]
    Y = CDS_grid[idx,:,:,sigma_fc_idx] / CDS_ref_grid[idx,:,:,sigma_fc_idx]
    
    N_eng_unique = N_eng_unique.astype(float)
    N_eng_unique[-1] -= 1e-6
    cs2 = ax_main.contour(
        X,
        Y,
        B[idx,:,:,sigma_fc_idx],
        levels=N_eng_unique,
        colors=colors[2],
        linestyles='dashed',
        linewidths=0.8,
        zorder=10,
        extend='max',
    )
    wing_frac_unique = wing_frac_unique.astype(float)
    wing_frac_unique[-1] -= 1e-6
    cs2 = ax_main.contour(
        X,
        Y,
        C[idx,:,:,sigma_fc_idx],
        levels=wing_frac_unique,
        colors=colors[2],
        linestyles='dashdot',
        linewidths=0.8,
        zorder=10,
        extend='max',
    )
    ax_main.contourf(
        X,
        Y,
        B[idx,:,:,sigma_fc_idx],
        levels=N_eng_unique,
        colors=colors[2],
        alpha=0.1,
        zorder=5,  # below your contour lines
    )

##

xlim = ax_main.get_xlim()
ylim = ax_main.get_ylim()

##

ax_main.plot([1, 2], [1, 2], color='black')

# ax_main.scatter(WMTO_grid[:,:,:,0] / WMTO_ref_grid[:,:,:,0], CDS_grid[:,:,:,0] / CDS_ref_grid[:,:,:,0], marker='.')

cvhull_points = np.array([
    WMTO_grid.flatten() / WMTO_ref_grid.flatten(),
    CDS_grid.flatten() / CDS_ref_grid.flatten(),
]).T
cvhull = ConvexHull(cvhull_points)
cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
polygon = Polygon(cvhull_vertices)
x_outline_polygon, y_outline_polygon = polygon.exterior.xy
ax_main.fill(
    x_outline_polygon, y_outline_polygon, color='grey', alpha=0.1, zorder=-1,
)

##

# cvhull_points = np.array([
#     WMTO_grid[:,:,:,0].flatten() / WMTO_ref_grid[:,:,:,0].flatten(),
#     CDS_grid[:,:,:,0].flatten() / CDS_ref_grid[:,:,:,0].flatten(),
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax_main.fill(
#     x_outline_polygon, y_outline_polygon, color='black', alpha=0.1, zorder=-1,
# )

##

# Create a proxy artist for the legend
space_patch = mpatches.Patch(color='grey', alpha=0.1, label='Space')

# Add it to your legend handles
legend_elements.append(space_patch)

##

# ax_main.scatter(WMTO_grid / WMTO_ref_grid, CDS_grid / CDS_ref_grid, marker='.')

# ax_main.set_xlabel("Relative weight (-)")
# ax_main.set_ylabel("Relative drag (-)")
ax_main.spines[['right', 'top']].set_visible(False)
ax_main.tick_params(axis='y', which='both', right=False, length=0)
ax_main.tick_params(axis='x', which='both', length=0)

ax_main.set_xlim(*xlim)
ax_main.set_ylim(*ylim)

ax_main.set_xlim(left=1)
ax_main.set_ylim(bottom=1)

ax_main.legend(handles=legend_elements, loc='upper left', frameon=False)

# ax_main.set_xscale('log')
# ax_main.set_yscale('log')

ax_main.set_aspect('equal')

#####

## Vary fcs_loc_unique and N_eng_unique

idx = 0

X = WMTO_grid[:,:,idx,sigma_fc_idx] / WMTO_ref_grid[:,:,idx,sigma_fc_idx]
Y = CDS_grid[:,:,idx,sigma_fc_idx] / CDS_ref_grid[:,:,idx,sigma_fc_idx]

fcs_loc_unique = fcs_loc_unique.astype(float)
fcs_loc_unique[-1] -= 1e-6
cs = ax_blue.contour(
    X,
    Y,
    A[:,:,idx,sigma_fc_idx],
    levels=fcs_loc_unique,
    colors=colors[0],
    linestyles='solid',
    linewidths=0.8,
    zorder=10,
    extend='max',
)
ax_blue.clabel(cs, fmt='%0.3f')
N_eng_unique = N_eng_unique.astype(float)
N_eng_unique[-1] -= 1e-6
cs2 = ax_blue.contour(
    X,
    Y,
    B[:,:,idx,sigma_fc_idx],
    levels=N_eng_unique,
    colors=colors[0],
    linestyles='dashed',
    linewidths=0.8,
    zorder=10,
    extend='max',
)
ax_blue.clabel(cs2, fmt='%0.3f')
ax_blue.contourf(
    X,
    Y,
    A[:,:,idx,sigma_fc_idx],
    levels=fcs_loc_unique,
    colors=colors[0],
    alpha=0.1,
    zorder=5,  # below your contour lines
)

# ax_blue.set_xlabel("Relative weight (-)")
# ax_blue.set_ylabel("Relative drag (-)")
ax_blue.spines[['right', 'top']].set_visible(False)
ax_blue.tick_params(axis='y', which='both', right=False, length=0)
ax_blue.tick_params(axis='x', which='both', length=0)
# ax_blue.set_aspect('equal')
add_margin(ax_blue, 0.05)

#####

## Vary fcs_loc_unique and wing_frac_idx

idx = 0

X = WMTO_grid[:,idx,:,sigma_fc_idx] / WMTO_ref_grid[:,idx,:,sigma_fc_idx]
Y = CDS_grid[:,idx,:,sigma_fc_idx] / CDS_ref_grid[:,idx,:,sigma_fc_idx]

fcs_loc_unique = fcs_loc_unique.astype(float)
fcs_loc_unique[-1] -= 1e-6
cs = ax_green.contour(
    X,
    Y,
    A[:,idx,:,sigma_fc_idx],
    levels=fcs_loc_unique,
    colors=colors[1],
    linestyles='solid',
    linewidths=0.8,
    zorder=10,
    extend='max',
)
ax_green.clabel(cs, fmt='%0.3f')
wing_frac_unique = wing_frac_unique.astype(float)
wing_frac_unique[-1] -= 1e-6
cs2 = ax_green.contour(
    X,
    Y,
    C[:,idx,:,sigma_fc_idx],
    levels=wing_frac_unique,
    colors=colors[1],
    linestyles='dashdot',
    linewidths=0.8,
    zorder=10,
    extend='max',
)
ax_green.clabel(cs2, fmt='%0.3f')
# Fill the area covered by the first contour set with red, alpha=0.5
ax_green.contourf(
    X,
    Y,
    A[:,idx,:,sigma_fc_idx],
    levels=fcs_loc_unique,
    colors=colors[1],
    alpha=0.1,
    zorder=5,  # below your contour lines
)

# ax_green.set_xlabel("Relative weight (-)")
# ax_green.set_ylabel("Relative drag (-)")
ax_green.spines[['right', 'top']].set_visible(False)
ax_green.tick_params(axis='y', which='both', right=False, length=0)
ax_green.tick_params(axis='x', which='both', length=0)
# ax_green.set_aspect('equal')
add_margin(ax_green, 0.05)

#####

## Vary N_eng_idx and wing_frac_idx

idx = -1

X = WMTO_grid[idx,:,:,sigma_fc_idx] / WMTO_ref_grid[idx,:,:,sigma_fc_idx]
Y = CDS_grid[idx,:,:,sigma_fc_idx] / CDS_ref_grid[idx,:,:,sigma_fc_idx]

N_eng_unique = N_eng_unique.astype(float)
N_eng_unique[-1] -= 1e-6
cs2 = ax_purple.contour(
    X,
    Y,
    B[idx,:,:,sigma_fc_idx],
    levels=N_eng_unique,
    colors=colors[2],
    linestyles='dashed',
    linewidths=0.8,
    zorder=10,
    extend='max',
)
ax_purple.clabel(cs2, fmt='%0.3f')
wing_frac_unique = wing_frac_unique.astype(float)
wing_frac_unique[-1] -= 1e-6
cs2 = ax_purple.contour(
    X,
    Y,
    C[idx,:,:,sigma_fc_idx],
    levels=wing_frac_unique,
    colors=colors[2],
    linestyles='dashdot',
    linewidths=0.8,
    zorder=10,
    extend='max',
)
ax_purple.clabel(cs2, fmt='%0.3f')
ax_purple.contourf(
    X,
    Y,
    B[idx,:,:,sigma_fc_idx],
    levels=N_eng_unique,
    colors=colors[2],
    alpha=0.1,
    zorder=5,  # below your contour lines
)

# ax_purple.set_xlabel("Relative weight (-)")
# ax_purple.set_ylabel("Relative drag (-)")
ax_purple.spines[['right', 'top']].set_visible(False)
ax_purple.tick_params(axis='y', which='both', right=False, length=0)
ax_purple.tick_params(axis='x', which='both', length=0)
# ax_purple.set_aspect('equal')
add_margin(ax_purple, 0.05)

#####

fig.text(0.5, 0.06, 'Relative weight (-)', ha='center', va='center')
fig.text(0.06, 0.5, 'Relative drag (-)', ha='center', va='center', rotation='vertical')

# plt.tight_layout()
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



