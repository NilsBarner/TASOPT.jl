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

xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'volume_results_narrowbody.xlsx'))

sheet_dict = {}
for sheet_name in xlsx.sheet_names:
    sheet_dict[sheet_name] = pd.read_excel(xlsx, sheet_name=sheet_name)

#%%

df = sheet_dict['results']

radius_unique = np.sort(df["radius"].unique())
AR_unique = np.sort(df["AR"].unique())
N_eng_unique = np.sort(df["N_eng"].unique())
Vspec_unique = np.sort(df["Vspec"].unique())

# Create 2D coordinate grids
A, B, C, D = np.meshgrid(radius_unique, AR_unique, N_eng_unique, Vspec_unique, indexing="ij")

# Create shape tuple for reshaping
shape = (
    len(radius_unique),
    len(AR_unique),
    len(N_eng_unique),
    len(Vspec_unique),
)

# Create index mapping for each dimension
radius_idx = {v: i for i, v in enumerate(radius_unique)}
AR_idx = {v: i for i, v in enumerate(AR_unique)}
N_eng_idx = {v: i for i, v in enumerate(N_eng_unique)}
Vspec_idx = {v: i for i, v in enumerate(Vspec_unique)}

# Initialize empty arrays
CDS_grid = np.full(shape, np.nan)
CDS_ref_grid = np.full(shape, np.nan)
WMTO_grid = np.full(shape, np.nan)
WMTO_ref_grid = np.full(shape, np.nan)

# Populate arrays
for _, row in df.iterrows():
    i = radius_idx[row["radius"]]
    j = AR_idx[row["AR"]]
    k = N_eng_idx[row["N_eng"]]
    l = Vspec_idx[row["Vspec"]]

    CDS_grid[i, j, k, l] = row["CDS"]
    CDS_ref_grid[i, j, k, l] = row["CDS_ref"]
    WMTO_grid[i, j, k, l] = row["WMTO"]
    WMTO_ref_grid[i, j, k, l] = row["WMTO_ref"]

#%%

def add_margin(ax, m=0.05):
    for a, s in [(ax.get_xlim, ax.set_xlim), (ax.get_ylim, ax.set_ylim)]:
        lo, hi = a(); r = hi - lo; s(lo - m*r, hi + m*r)

plt.rcParams['axes.prop_cycle'] = cycler('color', ['#206095', '#a8bd3a', '#871a5b', '#f66068', '#05341A', '#27a0cc', '#003c57', '#22d0b6', '#746cb1', '#A09FA0'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

legend_elements = [
    Line2D([0], [0], color='black', marker='o', label=r'Baseline $\left( N_\mathrm{eng}=2, \%_\mathrm{wing}=1, \sigma_\mathrm{FC}=2~\mathrm{kW/kg} \right)$'),
    Line2D([0], [0], color=colors[0], label=r'$\hat{x} \rightarrow \left[ 0,1 \right]$'),
    Line2D([0], [0], color=colors[1], label=r'$N_\mathrm{eng} \rightarrow \left[ 2,8 \right]$'),
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
    
Vspec_idx = 0

fig, ax = plt.subplots(figsize=(8, 6))

# Explicit axis limits needed when not using ax.scatter() due to nan's in some arrays
X = WMTO_grid / WMTO_ref_grid
Y = CDS_grid / CDS_ref_grid
ax.set_xlim(np.nanmin(X), np.nanmax(X))
ax.set_ylim(np.nanmin(Y), np.nanmax(Y))

xlim = ax.get_xlim()
ylim = ax.get_ylim()

##

# cvhull_points = np.array([
#     WMTO_grid.flatten() / WMTO_ref_grid.flatten(),
#     CDS_grid.flatten() / CDS_ref_grid.flatten(),
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax.fill(
#     x_outline_polygon, y_outline_polygon, color='grey', alpha=0.1, zorder=-1,
# )

##

# cvhull_points = np.array([
#     WMTO_grid[:,:,:,0].flatten() / WMTO_ref_grid[:,:,:,0].flatten(),
#     CDS_grid[:,:,:,0].flatten() / CDS_ref_grid[:,:,:,0].flatten(),
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax.fill(
#     x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
# )

##

# cvhull_points = np.array([
#     WMTO_grid[:,:,:,-1].flatten() / WMTO_ref_grid[:,:,:,-1].flatten(),
#     CDS_grid[:,:,:,-1].flatten() / CDS_ref_grid[:,:,:,-1].flatten(),
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax.fill(
#     x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
# )

##

# ax.scatter(WMTO_grid / WMTO_ref_grid, CDS_grid / CDS_ref_grid, marker='.')

# ax.scatter(WMTO_grid[2,:,:,:] / WMTO_ref_grid[2,:,:,:], CDS_grid[2,:,:,:] / CDS_ref_grid[2,:,:,:], marker='.', color='blue')
# ax.scatter(WMTO_grid[:,2,:,:] / WMTO_ref_grid[:,2,:,:], CDS_grid[:,2,:,:] / CDS_ref_grid[:,2,:,:], marker='.', color='green')
# ax.scatter(WMTO_grid[:,:,2,:] / WMTO_ref_grid[:,:,2,:], CDS_grid[:,:,2,:] / CDS_ref_grid[:,:,2,:], marker='.', color='red')
ax.scatter(WMTO_grid[:,:,:,2] / WMTO_ref_grid[:,:,:,2], CDS_grid[:,:,:,2] / CDS_ref_grid[:,:,:,2], marker='.', color='orange')

ax.set_xlabel("Relative weight (-)", labelpad=10)
ax.set_ylabel("Relative drag (-)", labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

# ax.set_xlim(*xlim)
# ax.set_ylim(*ylim)

# ax.set_xlim(left=1)
# ax.set_ylim(bottom=1)

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

idx_radius_des = -1  # inactive
idx_AR_des = 0
idx_N_eng_des = -1
idx_Vspec_des = 0

ax.scatter(
    WMTO_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, idx_Vspec_des] / WMTO_ref_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, idx_Vspec_des],
    CDS_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, idx_Vspec_des] / CDS_ref_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, idx_Vspec_des],
    marker='o',
    color='black',
    zorder=20,
)

for idx_radius in range(len(radius_unique)):
    # ax.plot(
    #     WMTO_grid[idx_radius, idx_AR_des, [0, -1], idx_Vspec_des] / WMTO_ref_grid[idx_radius, idx_AR_des, [0, -1], idx_Vspec_des],
    #     CDS_grid[idx_radius, idx_AR_des, [0, -1], idx_Vspec_des] / CDS_ref_grid[idx_radius, idx_AR_des, [0, -1], idx_Vspec_des],
    #     color=colors[0],
    #     marker='.',
    # )
    ax.plot(
        WMTO_grid[idx_radius, idx_AR_des, :, idx_Vspec_des] / WMTO_ref_grid[idx_radius, idx_AR_des, :, idx_Vspec_des],
        CDS_grid[idx_radius, idx_AR_des, :, idx_Vspec_des] / CDS_ref_grid[idx_radius, idx_AR_des, :, idx_Vspec_des],
        color=colors[0],
    )
    
    ax.scatter(
        WMTO_grid[idx_radius, idx_AR_des, :, idx_Vspec_des] / WMTO_ref_grid[idx_radius, idx_AR_des, :, idx_Vspec_des],
        CDS_grid[idx_radius, idx_AR_des, :, idx_Vspec_des] / CDS_ref_grid[idx_radius, idx_AR_des, :, idx_Vspec_des],
        color=colors[2],
        marker='.',
        zorder=10,
    )
    
# for idx_AR in range(len(AR_unique)):
ax.plot(
    WMTO_grid[idx_radius_des, :, idx_N_eng_des, idx_Vspec_des] / WMTO_ref_grid[idx_radius_des, :, idx_N_eng_des, idx_Vspec_des],
    CDS_grid[idx_radius_des, :, idx_N_eng_des, idx_Vspec_des] / CDS_ref_grid[idx_radius_des, :, idx_N_eng_des, idx_Vspec_des],
    color=colors[1],
    marker='.',
)
    
# # for idx_N_eng in range(len(N_eng_unique)):
# ax.scatter(
#     WMTO_grid[idx_radius_des, idx_AR_des, :, idx_Vspec_des] / WMTO_ref_grid[idx_radius_des, idx_AR_des, :, idx_Vspec_des],
#     CDS_grid[idx_radius_des, idx_AR_des, :, idx_Vspec_des] / CDS_ref_grid[idx_radius_des, idx_AR_des, :, idx_Vspec_des],
#     color=colors[2],
#     marker='.',
# )
    
# for idx_Vspec in range(len(Vspec_unique)):
ax.plot(
    WMTO_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, :] / WMTO_ref_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, :],
    CDS_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, :] / CDS_ref_grid[idx_radius_des, idx_AR_des, idx_N_eng_des, :],
    color=colors[3],
    marker='.',
)

# cvhull_points = np.array([
#     WMTO_grid[:,:,:,0].flatten() / WMTO_ref_grid[:,:,:,0].flatten(),
#     CDS_grid[:,:,:,0].flatten() / CDS_ref_grid[:,:,:,0].flatten(),
# ]).T
# cvhull = ConvexHull(cvhull_points)
# cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
# polygon = Polygon(cvhull_vertices)
# x_outline_polygon, y_outline_polygon = polygon.exterior.xy
# ax.fill(
#     x_outline_polygon, y_outline_polygon, color='black', alpha=0.2, zorder=-1,
# )

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



