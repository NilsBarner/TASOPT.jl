import os
import re
import sys
import ast
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

# df = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'data', 'volume_results_narrowbody_newesttttt.csv'))
# df_ref = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'data', 'volume_results_narrowbody_ref_newest.csv'))

# df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_101125.csv'))
# df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_101125.csv'))
df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_191125.csv'))
df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_191125.csv'))

#%%

index_py = np.vstack(df["index"].apply(ast.literal_eval).to_numpy()) - 1
df["index"] = [tuple(i) for i in index_py]
fcs_loc_unique = df["fcs_loc"].unique()
radius_unique = np.sort(df["radius"].unique())
AR_unique = np.sort(df["AR"].unique())
tdivc_scale_unique = np.sort(df["tdivc_scale"].unique())
N_eng_unique = np.sort(df["N_eng"].unique())
HTR_f_unique = np.sort(df["HTR_f"].unique())
# Vspec_unique = np.sort(df["Vspec"].unique())
l_fcs_unique = np.sort(df["l_fcs"].unique())
theta_floor_unique = np.sort(df["theta_floor"].unique())
fcs_fuselage_location_unique = np.sort(df["fcs_fuselage_location"].unique())
has_strut_unique = np.sort(df["has_strut"].unique())

#%%

fcs_loc_ref = df_ref["fcs_loc"][0]
radius_ref = df_ref["radius"][0]
AR_ref = df_ref["AR"][0]
tdivc_scale_ref = df_ref["tdivc_scale"][0]
N_eng_ref = df_ref["N_eng"][0]
HTR_f_ref = df_ref["HTR_f"][0]
l_fcs_ref = df_ref["l_fcs"][0]
theta_floor_ref = df_ref["theta_floor"][0]
fcs_fuselage_location_ref = df_ref["fcs_fuselage_location"][0]
has_strut_ref = df_ref["has_strut"][0]

CDS_ref = df_ref["CDS"]
WMTO_ref = df_ref["WMTO"]
Vol_wing_ref = df_ref["Vol_wing"]
seats_abreast_ref = df_ref["seats_abreast"]
PFEI_ref = df_ref["PFEI"]
Vol_nacelle_ref = df_ref["Vol_nacelle"]
y_centroid_wing_ref = df_ref["y_centroid_wing"]
L_fuse_ref = df_ref["L_fuse"]
y_centroid_nacelles_ref = df_ref["y_centroid_nacelles"]
span_ref = df_ref["span"]
x_centroid_fuse_ref = df_ref["x_centroid_fuse"]
Vol_fuse_ref = df_ref["Vol_fuse"]

radius_ref_idx = np.nanargmin(abs(radius_unique - radius_ref))
AR_ref_idx = np.nanargmin(abs(AR_unique - AR_ref))
tdivc_scale_ref_idx = np.nanargmin(abs(tdivc_scale_unique - tdivc_scale_ref))
N_eng_ref_idx = np.nanargmin(abs(N_eng_unique - N_eng_ref))
HTR_f_ref_idx = np.nanargmin(abs(HTR_f_unique - HTR_f_ref))
l_fcs_ref_idx = np.nanargmin(abs(l_fcs_unique - l_fcs_ref))
theta_floor_ref_idx = np.nanargmin(abs(theta_floor_unique - theta_floor_ref))
fcs_fuselage_location_ref_idx = 1
has_strut_ref_idx = 0

# sys.exit()

#%%

# Nacelles
mask_nacelle = df["index"].apply(lambda x: x[0] == 0)
df_filtered_nacelle = df[mask_nacelle]
print("df_filtered_nacelle =", df_filtered_nacelle)

# Wing
mask_wing = df["index"].apply(lambda x: x[0] == 1)
df_filtered_wing = df[mask_wing]

# Fuselage

mask_fuselage = df["index"].apply(lambda x: x[0] == 2)
df_filtered_fuselage = df[mask_fuselage]

# mask_fuselage_rear = df["index"].apply(
#     lambda x:
#     (x[0] == 2) and
#     # (x[2] == AR_ref_idx) and
#     (x[2] == 8) and
#     (x[3] == tdivc_scale_ref_idx) and
#     (x[4] == N_eng_ref_idx) and
#     (x[5] == HTR_f_ref_idx) and
#     (x[8] == 1) and
#     (x[9] == 0)
# )
# df_filtered_fuselage_rear = df[mask_fuselage_rear]

# mask_fuselage_underfloor = df["index"].apply(
#     lambda x:
#     (x[0] == 2) and
#     # (x[2] == AR_ref_idx) and
#     (x[2] == 8) and
#     (x[3] == tdivc_scale_ref_idx) and
#     (x[4] == N_eng_ref_idx) and
#     (x[5] == HTR_f_ref_idx) and
#     (x[8] == 0) and
#     (x[9] == 0)
# )
# df_filtered_fuselage_underfloor = df[mask_fuselage_underfloor]

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(df_filtered_nacelle['WMTO'], df_filtered_nacelle['CDS'], zorder=100)
ax.scatter(df_filtered_wing['WMTO'], df_filtered_wing['CDS'], zorder=100)
ax.scatter(df_filtered_fuselage['WMTO'], df_filtered_fuselage['CDS'], zorder=100)

plt.show()

'''
def closest_index_along_axis(arr, ref, axis):
    """Return the global index along `axis` of the value closest to `ref`."""
    # Collapse all other axes so we only search along the one axis
    # Use nanmean to handle potential NaNs gracefully
    arr_1d = np.nanmean(arr, axis=tuple(i for i in range(arr.ndim) if i != axis))
    diff = np.abs(arr_1d - ref)
    return int(np.nanargmin(diff))

# fcs_loc_ref_idx = closest_index_along_axis(fcs_loc_grid, fcs_loc_ref.to_numpy(), axis=0)
radius_ref_idx = closest_index_along_axis(radius_grid, radius_ref.to_numpy(), axis=1)
AR_ref_idx = closest_index_along_axis(AR_grid, AR_ref.to_numpy(), axis=2)
tdivc_scale_ref_idx = closest_index_along_axis(tdivc_scale_grid, tdivc_scale_ref.to_numpy(), axis=3)
N_eng_ref_idx = closest_index_along_axis(N_eng_grid, N_eng_ref.to_numpy(), axis=4)
HTR_f_ref_idx = closest_index_along_axis(HTR_f_grid, HTR_f_ref.to_numpy(), axis=5)
Vspec_ref_idx = closest_index_along_axis(Vspec_grid, Vspec_ref.to_numpy(), axis=6)
# fcs_fuselage_location_ref_idx = closest_index_along_axis(fcs_fuselage_location_grid, fcs_fuselage_location_ref.to_numpy(), axis=0)

radius_ref_closest = radius_unique[radius_ref_idx]
AR_ref_closest = AR_unique[AR_ref_idx]
tdivc_scale_ref_closest = tdivc_scale_unique[tdivc_scale_ref_idx]
N_eng_ref_closest = N_eng_unique[N_eng_ref_idx]
HTR_f_ref_closest = HTR_f_unique[HTR_f_ref_idx]
Vspec_ref_closest = Vspec_unique[Vspec_ref_idx]

print('radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, Vspec_ref_idx =', radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, Vspec_ref_idx)
print('radius_ref_closest, AR_ref_closest, tdivc_scale_ref_closest, N_eng_ref_closest, HTR_f_ref_closest, Vspec_ref_closest =', radius_ref_closest, AR_ref_closest, tdivc_scale_ref_closest, N_eng_ref_closest, HTR_f_ref_closest, Vspec_ref_closest)

idxs_nacelle_baseline = (0, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, -1, 0)
# idxs_wing_baseline = (1, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, -1, 0)
# idxs_fuselage_1_baseline = (2, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, -1, 0)
# idxs_fuselage_2_baseline = (2, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, -1, 1)
# mask = 

# fig, ax = plt.subplots(figsize=(8, 6))

# for i in range(len(radius_unique)):

#     idxs_nacelle_baseline = (slice(None), i, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, -1, -1)
#     ax.scatter(WMTO_grid[idxs_nacelle_baseline], CDS_grid[idxs_nacelle_baseline])
#     print(i, WMTO_grid[idxs_nacelle_baseline], CDS_grid[idxs_nacelle_baseline])
#     print()

# # ax.scatter(WMTO_ref, CDS_ref)
# # ax.scatter(WMTO_grid, CDS_grid)
# # ax.scatter(WMTO_grid[idxs_nacelle_baseline], CDS_grid[idxs_nacelle_baseline])
# # ax.scatter(WMTO_grid[idxs_wing_baseline], CDS_grid[idxs_wing_baseline])
# # ax.scatter(WMTO_grid[idxs_fuselage_1_baseline], CDS_grid[idxs_fuselage_1_baseline])
# # ax.scatter(WMTO_grid[idxs_fuselage_2_baseline], CDS_grid[idxs_fuselage_2_baseline])

# plt.show()

sys.exit()

#%%

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(WMTO_ref, CDS_ref)
ax.scatter(WMTO_grid[0, slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)], CDS_grid[0, slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)])
ax.scatter(WMTO_grid[1, slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)], CDS_grid[1, slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)])
ax.scatter(WMTO_grid[2, slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)], CDS_grid[2, slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)])

plt.show()
'''

sys.exit()

#%%



# fig, ax = plt.subplots(figsize=(8, 6))

# ax.scatter(WMTO_ref, CDS_ref, zorder=100)
# ax.scatter(df["WMTO"], df["CDS"])

# plt.show()

# x = 

minus1 = len(df.columns) - 1

idx_baseline_nacelle = (0, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, l_fcs_ref_idx, theta_floor_ref_idx, 1, 0)
idx_baseline_wing = (1, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, l_fcs_ref_idx, theta_floor_ref_idx, 1, 0)
idx_baseline_fuselage_1 = (2, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, l_fcs_ref_idx, theta_floor_ref_idx, 0, 0)
idx_baseline_fuselage_2 = (2, radius_ref_idx, AR_ref_idx, tdivc_scale_ref_idx, N_eng_ref_idx, HTR_f_ref_idx, l_fcs_ref_idx, theta_floor_ref_idx, 1, 0)

# # mask = df["index"].apply(lambda x: x[0] == 2)
# mask_nacelle = df["index"].apply(lambda x: x == idx_baseline_nacelle)
# mask_wing = df["index"].apply(lambda x: x == idx_baseline_wing)
# mask_fuselage_1 = df["index"].apply(lambda x: x == idx_baseline_fuselage_1)
# mask_fuselage_2 = df["index"].apply(lambda x: x == idx_baseline_fuselage_2)

# df_filtered_nacelle = df[mask_nacelle]
# df_filtered_wing = df[mask_wing]
# df_filtered_fuselage_1 = df[mask_fuselage_1]
# df_filtered_fuselage_2 = df[mask_fuselage_2]

# mask = df["index"].apply(lambda x: x[0] == 1).apply(lambda x: x[1] == radius_ref_idx)
# mask = df["index"].apply(lambda x: x[0] == 1).apply(lambda x: x[1] == radius_ref_idx)

# Nacelles

mask_HTR_f = df["index"].apply(
    lambda x:
    (x[0] == 0) and
    (x[1] == radius_ref_idx) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_HTR_f = df[mask_HTR_f]

mask_N_eng = df["index"].apply(
    lambda x:
    (x[0] == 0) and
    (x[1] == radius_ref_idx) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_N_eng = df[mask_N_eng]

# Wing

mask_AR = df["index"].apply(
    lambda x:
    (x[0] == 1) and
    (x[1] == radius_ref_idx) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_AR = df[mask_AR]

mask_tdivc_scale = df["index"].apply(
    lambda x:
    (x[0] == 1) and
    (x[1] == radius_ref_idx) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_tdivc_scale = df[mask_tdivc_scale]

# Fuselage

mask_radius_0 = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 0) and
    (x[9] == 0)
)
df_filtered_radius_0 = df[mask_radius_0]

mask_radius_1 = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_radius_1 = df[mask_radius_1]

mask_l_fcs = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    (x[1] == radius_ref_idx) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_l_fcs = df[mask_l_fcs]
print("df_filtered_l_fcs =", df_filtered_l_fcs)

mask_theta_floor = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    (x[1] == radius_ref_idx) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[8] == 0) and
    (x[9] == 0)
)
df_filtered_theta_floor = df[mask_theta_floor]

#%%

# --- Plot
fig, ax = plt.subplots(figsize=(8, 6))

ax.set_aspect('equal')

ax.scatter(1, 1, color="black", marker='o', label="Reference", zorder=1000)

ax.plot(df_filtered_N_eng["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_N_eng["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[0], marker='.', linestyle='solid', zorder=100)
ax.plot(df_filtered_HTR_f["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_HTR_f["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[0], marker='.', linestyle='dashed', zorder=100)

ax.plot(df_filtered_AR["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_AR["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[1], marker='.', linestyle='solid')
ax.plot(df_filtered_tdivc_scale["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_tdivc_scale["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[1], marker='.', linestyle='dashed')

# ax.plot(df_filtered_radius_0["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_radius_0["CDS"].to_numpy() / CDS_ref.to_numpy(), marker='.', color=colors[2], linestyle='solid')
ax.plot(df_filtered_radius_1["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_radius_1["CDS"].to_numpy() / CDS_ref.to_numpy(), marker='.', color=colors[2], linestyle='solid')
ax.plot(df_filtered_l_fcs["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_l_fcs["CDS"].to_numpy() / CDS_ref.to_numpy(), marker='.', color=colors[2], linestyle='dashed')
ax.plot(df_filtered_theta_floor["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_theta_floor["CDS"].to_numpy() / CDS_ref.to_numpy(), marker='s', color=colors[2], linestyle='dashdot')

# ax.scatter(df_filtered_nacelle["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_nacelle["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[0], marker='.', linestyle='solid')
# ax.scatter(df_filtered_wing["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_wing["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[1], marker='.', linestyle='solid')
# ax.scatter(df_filtered_fuselage["WMTO"].to_numpy() / WMTO_ref.to_numpy(), df_filtered_fuselage["CDS"].to_numpy() / CDS_ref.to_numpy(), color=colors[2], marker='.', linestyle='solid')

# =============================================================================
import alphashape
from alpha_shape_creator import ensure_xy_array, add_shapely_polygon

points_2d = np.array([
    df_filtered_nacelle["WMTO"].to_numpy().flatten() / WMTO_ref.to_numpy().flatten(),
    df_filtered_nacelle["CDS"].to_numpy().flatten() / CDS_ref.to_numpy().flatten(),
]).T
pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
alpha_shape = alphashape.alphashape(pts, 1)
add_shapely_polygon(
    ax, alpha_shape, facecolor=colors[0], edgecolor='None', alpha=0.2,
)

points_2d = np.array([
    df_filtered_wing["WMTO"].to_numpy().flatten() / WMTO_ref.to_numpy().flatten(),
    df_filtered_wing["CDS"].to_numpy().flatten() / CDS_ref.to_numpy().flatten(),
]).T
pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
alpha_shape = alphashape.alphashape(pts, 20)
add_shapely_polygon(
    ax, alpha_shape, facecolor=colors[1], edgecolor='None', alpha=0.2,
)

points_2d = np.array([
    df_filtered_fuselage["WMTO"].to_numpy().flatten() / WMTO_ref.to_numpy().flatten(),
    df_filtered_fuselage["CDS"].to_numpy().flatten() / CDS_ref.to_numpy().flatten(),
]).T
pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
alpha_shape = alphashape.alphashape(pts, 10)
add_shapely_polygon(
    ax, alpha_shape, facecolor=colors[2], edgecolor='None', alpha=0.2,
)
# =============================================================================

ax.set_xlim(0.75, 1.4)
ax.set_ylim(0.875, 1.6)
# ax.set_xlim(right=1.41)
# ax.set_ylim(top=1.6)
ax.plot([0.9, 1.5], [0.9, 1.5], color='black', alpha=0.2)

ax.set_xlabel('Relative weight (-)', labelpad=10)
ax.set_ylabel('Relative drag (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

# build the handle lists exactly as before
handles_nacelle = [
    Line2D([0], [0], color=colors[0], linestyle='solid', marker='.', alpha=1.0, label=r'$N_\mathrm{eng}$'),
    Line2D([0], [0], color=colors[0], linestyle='dashed', marker='.', alpha=1.0, label=r'$\mathrm{HTR}_\mathrm{fan}$'),
]
handles_wing = [
    Line2D([0], [0], color=colors[1], linestyle='solid', marker='.', alpha=1.0, label=r'AR'),
    Line2D([0], [0], color=colors[1], linestyle='dashed', marker='.', alpha=1.0, label=r't/c'),
]
handles_fuselage = [
    Line2D([0], [0], color=colors[2], linestyle='solid', marker='.', alpha=1.0, label=r'$R_\mathrm{fuse}$'),
    Line2D([0], [0], color=colors[2], linestyle='dashed', marker='.', alpha=1.0, label=r'$L_\mathrm{FCS}$'),
    Line2D([0], [0], color=colors[2], linestyle='dashdot', marker='.', alpha=1.0, label=r'$\theta_\mathrm{floor}$'),
]

# Reserve some room on the right for the legends (tune this number if needed)
fig.subplots_adjust(right=0.50)

# Place three separate legends in figure coordinates (bbox_transform=fig.transFigure)
# y positions chosen to stack them vertically on the right
leg_nac = fig.legend(
    handles=handles_nacelle,
    title="Nacelle",
    loc='upper left',
    bbox_to_anchor=(0.75, 0.925),           # (x, y) in FIGURE fraction coords
    bbox_transform=fig.transFigure,
    frameon=False,
)
fig.add_artist(leg_nac)

leg_wing = fig.legend(
    handles=handles_wing,
    title="Wing",
    loc='upper left',
    bbox_to_anchor=(0.75, 0.675),
    bbox_transform=fig.transFigure,
    frameon=False,
)
fig.add_artist(leg_wing)

leg_fus = fig.legend(
    handles=handles_fuselage,
    title="Fuselage",
    loc='upper left',
    bbox_to_anchor=(0.75, 0.425),
    bbox_transform=fig.transFigure,
    frameon=False,
)
fig.add_artist(leg_fus)
# =============================================================================

# =============================================================================
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Create zoomed inset
axins = inset_axes(ax, width="115%", height="115%", loc='upper left',
                   bbox_to_anchor=(0.7625, 1.3, 0.3, 0.3),
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
axins.set_xlim(0.975, 1.05)
axins.set_ylim(0.975, 1.05)
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
# =============================================================================

plt.tight_layout()
plt.subplots_adjust(right=0.75)
# fig.tight_layout(rect=[0, 0, 0.75, 1])

# plt.savefig('weight_drag_volume.svg', format='svg')

plt.show()

#%%

print('Hello world')

fig, ax = plt.subplots(figsize=(8, 6))

# ax.scatter(df_filtered_fuselage_rear["x_centroid_fuse"]/df_filtered_fuselage_rear["length"], df_filtered_fuselage_rear["Vol_fuse"], c=df_filtered_fuselage_rear["seats_abreast"])
# ax.scatter(df_filtered_fuselage_underfloor["x_centroid_fuse"]/df_filtered_fuselage_underfloor["length"], df_filtered_fuselage_underfloor["Vol_fuse"], c=df_filtered_fuselage_underfloor["seats_abreast"])
# # ax.tricontourf(df_filtered_fuselage_rear["x_centroid_fuse"], df_filtered_fuselage_rear["Vol_fuse"], df_filtered_fuselage_rear["seats_abreast"])
# # ax.tricontourf(df_filtered_fuselage_underfloor["x_centroid_fuse"], df_filtered_fuselage_underfloor["Vol_fuse"], df_filtered_fuselage_underfloor["seats_abreast"])

# ax.scatter(df_filtered_fuselage_rear["x_centroid_fuse"]/df_filtered_fuselage_rear["length"], df_filtered_fuselage_rear["Vol_fuse"], c=df_filtered_fuselage_rear["radius"])
# ax.scatter(df_filtered_fuselage_underfloor["x_centroid_fuse"]/df_filtered_fuselage_underfloor["length"], df_filtered_fuselage_underfloor["Vol_fuse"], c=df_filtered_fuselage_underfloor["radius"])
# # ax.tricontourf(df_filtered_fuselage_rear["x_centroid_fuse"], df_filtered_fuselage_rear["Vol_fuse"], df_filtered_fuselage_rear["radius"])
# # ax.tricontourf(df_filtered_fuselage_underfloor["x_centroid_fuse"], df_filtered_fuselage_underfloor["Vol_fuse"], df_filtered_fuselage_underfloor["radius"])

# ax.scatter(df_filtered_fuselage_rear["x_centroid_fuse"]/df_filtered_fuselage_rear["length"], df_filtered_fuselage_rear["Vol_fuse"], c=df_filtered_fuselage_rear["l_fcs"])
# ax.scatter(df_filtered_fuselage_underfloor["x_centroid_fuse"]/df_filtered_fuselage_underfloor["length"], df_filtered_fuselage_underfloor["Vol_fuse"], c=df_filtered_fuselage_underfloor["l_fcs"])
# # ax.tricontourf(df_filtered_fuselage_rear["x_centroid_fuse"], df_filtered_fuselage_rear["Vol_fuse"], df_filtered_fuselage_rear["l_fcs"])
# # ax.tricontourf(df_filtered_fuselage_underfloor["x_centroid_fuse"], df_filtered_fuselage_underfloor["Vol_fuse"], df_filtered_fuselage_underfloor["l_fcs"])

ax.scatter(df_filtered_fuselage_rear["x_centroid_fuse"]/df_filtered_fuselage_rear["length"], df_filtered_fuselage_rear["Vol_fuse"], c=df_filtered_fuselage_rear["theta_floor"])
ax.scatter(df_filtered_fuselage_underfloor["x_centroid_fuse"]/df_filtered_fuselage_underfloor["length"], df_filtered_fuselage_underfloor["Vol_fuse"], c=df_filtered_fuselage_underfloor["theta_floor"])
# ax.tricontourf(df_filtered_fuselage_rear["x_centroid_fuse"], df_filtered_fuselage_rear["Vol_fuse"], df_filtered_fuselage_rear["theta_floor"])
# ax.tricontourf(df_filtered_fuselage_underfloor["x_centroid_fuse"], df_filtered_fuselage_underfloor["Vol_fuse"], df_filtered_fuselage_underfloor["theta_floor"])
"""
# =============================================================================
import alphashape
from alpha_shape_creator import ensure_xy_array, add_shapely_polygon

points_2d = np.array([
    df_filtered_fuselage_rear["x_centroid_fuse"].to_numpy().flatten() / df_filtered_fuselage_rear["length"].to_numpy().flatten(),
    df_filtered_fuselage_rear["Vol_fuse"].to_numpy().flatten(),
]).T
pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
alpha_shape = alphashape.alphashape(pts, 0.001)
add_shapely_polygon(
    ax, alpha_shape, facecolor=colors[0], edgecolor='None', alpha=0.2,
)

points_2d = np.array([
    df_filtered_fuselage_underfloor["x_centroid_fuse"].to_numpy().flatten() / df_filtered_fuselage_underfloor["length"].to_numpy().flatten(),
    df_filtered_fuselage_underfloor["Vol_fuse"].to_numpy().flatten(),
]).T
pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
alpha_shape = alphashape.alphashape(pts, 0.001)
add_shapely_polygon(
    ax, alpha_shape, facecolor=colors[1], edgecolor='None', alpha=0.2,
)

seats_abreast_unique = np.unique(df["seats_abreast"].to_numpy())
for i, seats_abreast in enumerate(seats_abreast_unique):
    mask_seats_abreast = df["seats_abreast"].apply(lambda x: x == seats_abreast)
    df_filtered_seats_abreast = df[mask_seats_abreast]
    ax.scatter(
        df_filtered_seats_abreast["x_centroid_fuse"].to_numpy().flatten() / df_filtered_seats_abreast["length"].to_numpy().flatten(), df_filtered_seats_abreast["Vol_fuse"], color=colors[i]
    )
# =============================================================================

ax.set_xlabel('Location of centroid expressed as fraction of fuselage length (-)', labelpad=10)
ax.set_ylabel('Available / required FCS volume (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_elements = [
    Patch(facecolor=colors[0], edgecolor='None', alpha=0.2, label='Aft FCS'),
    Patch(facecolor=colors[1], edgecolor='None', alpha=0.2, label='Underfloor FCS'),
    Line2D([0], [0], color=colors[0], marker='o', alpha=1.0, linewidth=0, label=r'6 seats abreast'),
    Line2D([0], [0], color=colors[1], marker='o', alpha=1.0, linewidth=0, label=r'7 seats abreast'),
    Line2D([0], [0], color=colors[2], marker='o', alpha=1.0, linewidth=0, label=r'8 seats abreast'),
    Line2D([0], [0], color=colors[3], marker='o', alpha=1.0, linewidth=0, label=r'9 seats abreast'),
    Line2D([0], [0], color=colors[4], marker='o', alpha=1.0, linewidth=0, label=r'10 seats abreast'),
]
leg_main = ax.legend(handles=legend_elements, loc='upper right', frameon=False)
ax.add_artist(leg_main)
"""
plt.show()

sys.exit()

#%%

P_req = 20e6
Vspec_curr = 0.5e6
Vol_fuse_req = P_req / Vspec_curr
# P_tot = 20e6

fig, (axL, axR) = plt.subplots(1,2,sharey=True,figsize=(8,6),gridspec_kw={'wspace':0.075})
import alphashape
from alpha_shape_creator import ensure_xy_array, add_shapely_polygon
def _plot_region(ax,xmin,xmax):
    maskR = (df_filtered_fuselage_rear["x_centroid_fuse"]/df_filtered_fuselage_rear["length"]>=xmin) & (df_filtered_fuselage_rear["x_centroid_fuse"]/df_filtered_fuselage_rear["length"]<=xmax)
    pts = np.array([
        df_filtered_fuselage_rear["x_centroid_fuse"][maskR].to_numpy()/df_filtered_fuselage_rear["length"][maskR].to_numpy(),
        df_filtered_fuselage_rear["Vol_fuse"][maskR].to_numpy() / Vol_fuse_req
        # P_tot / df_filtered_fuselage_rear["Vol_fuse"][maskR].to_numpy() / 1e6
    ]).T
    if pts.size:
        add_shapely_polygon(
            ax, alphashape.alphashape(ensure_xy_array(pts),0.00001), facecolor=colors[0], edgecolor='None', alpha=0.3
        )
    maskU = (df_filtered_fuselage_underfloor["x_centroid_fuse"]/df_filtered_fuselage_underfloor["length"]>=xmin)&(df_filtered_fuselage_underfloor["x_centroid_fuse"]/df_filtered_fuselage_underfloor["length"]<=xmax)
    pts = np.array([
        df_filtered_fuselage_underfloor["x_centroid_fuse"][maskU].to_numpy()/df_filtered_fuselage_underfloor["length"][maskU].to_numpy(),
        df_filtered_fuselage_underfloor["Vol_fuse"][maskU].to_numpy() / Vol_fuse_req
        # P_tot / df_filtered_fuselage_underfloor["Vol_fuse"][maskU].to_numpy() / 1e6
    ]).T
    if pts.size:
        add_shapely_polygon(
            ax, alphashape.alphashape(ensure_xy_array(pts),0.00001), facecolor=colors[1], edgecolor='None', alpha=0.3
        )
    seats_abreast_unique = np.unique(df["seats_abreast"].to_numpy())
    for i, seats_abreast in enumerate(seats_abreast_unique):
        dfm = df[df["seats_abreast"]==seats_abreast]
        m = (dfm["x_centroid_fuse"]/dfm["length"]>=xmin)&(dfm["x_centroid_fuse"]/dfm["length"]<=xmax)
        if m.any():
            ax.scatter(
                (dfm["x_centroid_fuse"]/dfm["length"])[m].to_numpy(),
                dfm["Vol_fuse"][m].to_numpy() / Vol_fuse_req,
                # P_tot / dfm["Vol_fuse"][m].to_numpy() / 1e6,
                color=colors[i], marker='.',
            )
_plot_region(axL, 0.35, 0.425)
_plot_region(axR, 0.6, 0.65)
axL.set_xlim(0.35,0.41)
axR.set_xlim(0.605,0.655)
# axL.set_ylabel('Required FCS power density (kW/l)', labelpad=10)
axL.set_ylabel('Available / required FCS volume (-)', labelpad=10)
axL.spines[['right','top']].set_visible(False)
axR.spines[['right','top','left']].set_visible(False)
axL.tick_params(axis='y', which='both', right=False, length=0)
axL.tick_params(axis='x', which='both', length=0)
axR.tick_params(axis='y', which='both', right=False, length=0)
axR.tick_params(axis='x', which='both', length=0)
d = .015
kwargs = dict(transform=axL.transAxes, color='k', clip_on=False)
axL.plot((1-d,1+d),(-d,+d), **kwargs)
# axL.plot((1-d,1+d),(1-d,1+d), **kwargs)
kwargs = dict(transform=axR.transAxes, color='k', clip_on=False)
axR.plot((-d,+d),(-d,+d), **kwargs)
# axR.plot((-d,+d),(1-d,1+d), **kwargs)
legend_elements = [
    Patch(facecolor=colors[0], edgecolor='None', alpha=0.3, label='Aft FCS'),
    Patch(facecolor=colors[1], edgecolor='None', alpha=0.3, label='Underfloor FCS'),
    Line2D([0],[0], color=colors[0], marker='.', alpha=1.0, linewidth=0, label=r'6 seats abreast'),
    Line2D([0],[0], color=colors[1], marker='.', alpha=1.0, linewidth=0, label=r'7 seats abreast'),
    Line2D([0],[0], color=colors[2], marker='.', alpha=1.0, linewidth=0, label=r'8 seats abreast'),
    Line2D([0],[0], color=colors[3], marker='.', alpha=1.0, linewidth=0, label=r'9 seats abreast'),
    Line2D([0],[0], color=colors[4], marker='.', alpha=1.0, linewidth=0, label=r'10 seats abreast'),
]
leg = axR.legend(handles=legend_elements, loc='upper right', frameon=False)
leg.set_zorder(100)

fig.text(0.5, 0.0, 'Location of centroid expressed as fraction of fuselage length (-)', horizontalalignment='center')
# plt.savefig('fuselage_centroid_vs_volume.svg', format='svg')
plt.show()

#%%

# # =============================================================================
# for i, radius in enumerate(radius_unique):
    
#     mask_radius_curr = df["index"].apply(
#         lambda x:
#         (x[0] == 2) and
#         # (x[1] == i) and
#         (x[2] == AR_ref_idx) and
#         (x[3] == tdivc_scale_ref_idx) and
#         (x[4] == N_eng_ref_idx) and
#         (x[5] == HTR_f_ref_idx) and
#         (x[6] == l_fcs_ref_idx) and
#         (x[7] == theta_floor_ref_idx) and
#         (x[8] == 0) and
#         (x[9] == 0)
#     )
#     df_filtered_radius_curr = df[mask_radius_curr]
    
#     ax.plot(df_filtered_radius_curr["x_centroid_fuse"], df_filtered_radius_curr["Vol_fuse"])

# for i, l_fcs in enumerate(l_fcs_unique):
    
#     mask_l_fcs_curr = df["index"].apply(
#         lambda x:
#         (x[0] == 2) and
#         (x[1] == radius_ref) and
#         (x[2] == AR_ref_idx) and
#         (x[3] == tdivc_scale_ref_idx) and
#         (x[4] == N_eng_ref_idx) and
#         (x[5] == HTR_f_ref_idx) and
#         # (x[6] == l_fcs_ref_idx) and
#         (x[7] == theta_floor_ref_idx) and
#         (x[8] == 0) and
#         (x[9] == 0)
#     )
#     df_filtered_l_fcs_curr = df[mask_l_fcs_curr]
    
#     ax.plot(df_filtered_l_fcs_curr["x_centroid_fuse"], df_filtered_l_fcs_curr["Vol_fuse"])
# # =============================================================================


