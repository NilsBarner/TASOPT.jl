import os
import re
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
from cycler import cycler
from matplotlib import gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
from shapely import Polygon

from matplotlib_custom_settings import *

#%% Load all files of given name pattern

folder = os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'mass', 'data')  # or any other folder path

# Find all CSV files that start with 'combined_' and end with '_<number>.csv'
csv_files = glob.glob(os.path.join(folder, 'point_load_study_results_211125_[1-8]*.csv'))

# Sort numerically by the trailing number (so ..._2.csv comes before ..._10.csv)
def extract_number(filename):
    match = re.search(r'_(\d+)\.csv$', filename)
    return int(match.group(1)) if match else -1

csv_files = sorted(csv_files, key=extract_number)
sigma_fcs_list_sorted = [extract_number(csv_file) for csv_file in csv_files]

# Read all into a list of DataFrames
dfs = [pd.read_csv(f) for f in csv_files]

print(f"Loaded {len(dfs)} matching CSV files.")

#%% Reference and baseline aircraft

mask_bl = dfs[0]["nacelle_frac"].apply(lambda x: x == 1.0)
df_bl = dfs[0][mask_bl]
df_ref = pd.read_csv(os.path.join(folder, 'point_load_study_results_211125_ref.csv'))

#%% Plot bands of PFEI for nacelle, wing, and fuselage

from paretoset import paretoset

lw_plot = 2
alpha_pareto = 1
pareto_labels = ["Nacelle", "Wing", "Fuselage"]

fig, ax = plt.subplots(figsize=(10,8))

ax.axhline(1.0, color='k', linestyle='dashed', alpha=0.1, zorder=-100)

# Custom order to have nacelle -> wing -> fuselage (see mass_study_exec_wing.jl)
for i, df in enumerate([dfs[7], dfs[5], dfs[6]]):
    
    # Extract data for reusability
    sigma_fcs = df["sigma_fcs"].to_numpy() / 1e3
    PFEI_rel = df["PFEI"].to_numpy() / df_ref["PFEI"].to_numpy()
    
    # Nacelle data is a single line
    if i == 0:
        ax.plot(sigma_fcs, PFEI_rel, marker='.', color=colors[i], alpha=0.1)
    else:
        ax.scatter(sigma_fcs, PFEI_rel, marker='.', facecolor=colors[i], edgecolor='None', alpha=0.1)
    
    # Fill area between minima and maxima for each specific power
    
    # Extract minima and maxima for each specific power
    g = df.groupby("sigma_fcs")["PFEI"].agg(["min", "max"])
    # Convert to 1D numpy arrays
    sigma_vals = g.index.to_numpy() / 1e3
    PFEI_min = g["min"].to_numpy()
    PFEI_max = g["max"].to_numpy()
    ax.fill_between(x=sigma_vals, y1=PFEI_min / df_ref["PFEI"].to_numpy(), y2=PFEI_max / df_ref["PFEI"].to_numpy(), alpha=0.1, facecolor=colors[i], edgecolor='None', zorder=-i)
    
    # Plot pareto fronts
    points_2d = np.array([
        sigma_fcs,
        PFEI_rel,
    ]).T
    pareto_front = points_2d[paretoset(points_2d, sense=["min", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[i], linewidth=lw_plot, alpha=alpha_pareto, label=pareto_labels[i])

ax.set_xlabel('FCS specific power (kW/kg)', labelpad=10)
ax.set_ylabel('Relative PFEI (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
    
ax.legend(frameon=False, loc='upper right')
# plt.savefig("PFEI_vs_sigma.png", format='png', dpi=600)
# plt.savefig("PFEI_vs_sigma.svg", format='svg')

plt.show()

# sys.exit()

#%% Plot weight and drag for varying weight location and distribution, but fixed specific power

fig, ax = plt.subplots(figsize=(10,8))

# Plot lines for one-variation-at-a-time cases
ax.scatter(df_bl["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_bl["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='o', color='black', zorder=100)
ax.plot(dfs[0]["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), dfs[0]["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='.')
ax.plot(dfs[1]["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), dfs[1]["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='.')
ax.plot(dfs[2]["WMTO"].to_numpy()[:-1] / df_ref["WMTO"].to_numpy(), dfs[2]["CDS"].to_numpy()[:-1] / df_ref["CDS"].to_numpy(), marker='.')
ax.plot(dfs[3]["WMTO"].to_numpy()[:-1] / df_ref["WMTO"].to_numpy(), dfs[3]["CDS"].to_numpy()[:-1] / df_ref["CDS"].to_numpy(), marker='.')

# Plot grid for varying weight distribution

fcs_loc_unique = np.sort(dfs[4]["fcs_loc"].unique())
span_loc_unique = np.sort(dfs[4]["span_loc"].unique())
wing_frac_unique = np.sort(dfs[4]["wing_frac"].unique())

for fcs_loc in fcs_loc_unique:
    mask_temp = np.isclose(dfs[4]["fcs_loc"], fcs_loc)
    df_temp = dfs[4][mask_temp]
    
    ax.plot(df_temp["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_temp["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), color=colors[4], alpha=0.1, zorder=-20)
    
for span_loc in span_loc_unique:
    mask_temp = np.isclose(dfs[4]["span_loc"], span_loc)
    df_temp = dfs[4][mask_temp]
    
    ax.plot(df_temp["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_temp["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), color=colors[4], alpha=0.1, zorder=-20)

"""
###
# (x, y, z) = (
#     dfs[2]["WMTO"].to_numpy()[:-1] / df_ref["WMTO"].to_numpy(),
#     dfs[2]["CDS"].to_numpy()[:-1] / df_ref["CDS"].to_numpy(),
#     dfs[2]["fcs_loc"],
# )
(x, y, z) = (
    dfs[2]["WMTO"].to_numpy()[:-1],
    dfs[2]["CDS"].to_numpy()[:-1],
    dfs[2]["fcs_loc"],
)

for xi, yi, zi in zip(x, y, z):
    ax.text(
        xi, yi,
        f"{zi:.2g}",            # <-- one significant figure
        ha='left', va='top'  # tweak as desired
    )

# (x, y, z) = (
#     dfs[3]["WMTO"].to_numpy()[:-1] / df_ref["WMTO"].to_numpy(),
#     dfs[3]["CDS"].to_numpy()[:-1] / df_ref["CDS"].to_numpy(),
#     dfs[3]["span_loc"],
# )
(x, y, z) = (
    dfs[3]["WMTO"].to_numpy()[:-1],
    dfs[3]["CDS"].to_numpy()[:-1],
    dfs[3]["span_loc"],
)

for xi, yi, zi in zip(x, y, z):
    ax.text(
        xi, yi,
        f"{zi:.2g}",            # <-- one significant figure
        ha='left', va='top'  # tweak as desired
    )
###
"""

ax.set_xlabel('Relative weight (-)', labelpad=10)
ax.set_ylabel('Relative drag (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_handles = [
    Line2D([0], [0], linewidth=0.0, color='k', linestyle='solid', marker='o', alpha=1.0, label=r'Baseline'),
    Line2D([0], [0], color=colors[0], linestyle='solid', marker='.', alpha=1.0, label=r'Nacelle $\rightarrow$ fuselage'),
    Line2D([0], [0], color=colors[1], linestyle='solid', marker='.', alpha=1.0, label=r'Nacelle $\rightarrow$ wing'),
    Line2D([0], [0], color=colors[2], linestyle='solid', marker='.', alpha=1.0, label='Along fuselage'),
    Line2D([0], [0], color=colors[3], linestyle='solid', marker='.', alpha=1.0, label='Along wing'),
    Patch(facecolor='none', edgecolor=colors[4], alpha=0.1, label='Along wing and fuselage (50\% each)'),
]
legend = fig.legend(
    handles=legend_handles,
    loc='upper left',
    bbox_to_anchor=(0.13, 1),
    bbox_transform=fig.transFigure,
    frameon=False,
)
fig.add_artist(legend)
    
# plt.savefig("weight_vs_drag_mass_study_1to5.png", format='png', dpi=600)
# plt.savefig("weight_vs_drag_mass_study_1to5.svg", format='svg')
plt.show()

# sys.exit()

#%% Plot weight and drag for varying weight location and specific power

fig, ax = plt.subplots(figsize=(10,8))
    
# Plot line for nacelle
ax.scatter(df_bl["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_bl["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='o', color='black', zorder=100)
# ax.scatter(dfs[5]["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), dfs[5]["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='.')
# ax.scatter(dfs[6]["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), dfs[6]["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='.')
ax.plot(dfs[7]["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), dfs[7]["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), marker='.', color=colors[0])

# Plot grids for wing and fuselage

span_loc_unique = np.sort(dfs[5]["span_loc"].unique())
sigma_fcs_unique = np.sort(dfs[5]["sigma_fcs"].unique())

for span_loc in span_loc_unique:
    mask_temp = np.isclose(dfs[5]["span_loc"], span_loc)
    df_temp = dfs[5][mask_temp]
    
    ax.plot(df_temp["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_temp["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), color=colors[1], alpha=0.5, zorder=-10)
    
for sigma_fcs in sigma_fcs_unique:
    mask_temp = np.isclose(dfs[5]["sigma_fcs"], sigma_fcs)
    df_temp = dfs[5][mask_temp]
    
    ax.plot(df_temp["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_temp["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), color=colors[1], alpha=0.5, zorder=-10)
    

fcs_loc_unique = np.sort(dfs[6]["fcs_loc"].unique())
sigma_fcs_unique = np.sort(dfs[6]["sigma_fcs"].unique())
    
for fcs_loc in fcs_loc_unique:
    mask_temp = np.isclose(dfs[6]["fcs_loc"], fcs_loc)
    df_temp = dfs[6][mask_temp]
    
    ax.plot(df_temp["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_temp["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), color=colors[2], alpha=0.5, zorder=-20)
    
for sigma_fcs in sigma_fcs_unique:
    mask_temp = np.isclose(dfs[6]["sigma_fcs"], sigma_fcs)
    df_temp = dfs[6][mask_temp]
    
    ax.plot(df_temp["WMTO"].to_numpy() / df_ref["WMTO"].to_numpy(), df_temp["CDS"].to_numpy() / df_ref["CDS"].to_numpy(), color=colors[2], alpha=0.5, zorder=-20)


ax.set_xlabel('Relative weight (-)', labelpad=10)
ax.set_ylabel('Relative drag (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_handles = [
    Line2D([0], [0], linewidth=0.0, color='k', linestyle='solid', marker='o', alpha=1.0, label=r'Baseline'),
    Line2D([0], [0], color=colors[0], linestyle='solid', marker='.', alpha=1.0, label=r'Nacelle'),
    # Line2D([0], [0], color=colors[1], linestyle='solid', marker='.', alpha=1.0, label=r'Wing'),
    # Line2D([0], [0], color=colors[2], linestyle='solid', marker='.', alpha=1.0, label='Fuselage'),
    # Line2D([0], [0], color=colors[3], linestyle='solid', marker='.', alpha=1.0, label='Along wing'),
    Patch(facecolor='none', edgecolor=colors[1], alpha=0.5, label='Wing'),
    Patch(facecolor='none', edgecolor=colors[2], alpha=0.5, label='Fuselage'),
]
legend = fig.legend(
    handles=legend_handles,
    loc='upper left',
    bbox_to_anchor=(0.13, 1),
    bbox_transform=fig.transFigure,
    frameon=False,
)
fig.add_artist(legend)
    
# plt.savefig("weight_vs_drag_mass_study_6to8.png", format='png', dpi=600)
# plt.savefig("weight_vs_drag_mass_study_6to8.svg", format='svg')
plt.show()

sys.exit()

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



