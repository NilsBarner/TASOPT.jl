import os
import re
import sys
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import griddata
from cycler import cycler
from matplotlib import gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
from shapely import Polygon

from matplotlib_custom_settings import *

#%% Data import

df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_111125.csv'))
df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_111125.csv'))


#%%
# =============================================================================

nu_fcs_range = np.linspace(0.2, 5.5, 10) * 1e6

sigma_fcs_nacelle_list = []
sigma_fcs_wing_list = []
sigma_fcs_fuselage_rear_list = []
sigma_fcs_fuselage_underfloor_list = []

Vol_nacelle_list = []
Vol_wing_list = []
Vol_fuselage_rear_list = []
Vol_fuselage_underfloor_list = []

PFEI_nacelle_list = []
PFEI_wing_list = []
PFEI_fuselage_rear_list = []
PFEI_fuselage_underfloor_list = []

WMTO_nacelle_list = []
WMTO_wing_list = []
WMTO_fuselage_rear_list = []
WMTO_fuselage_underfloor_list = []

P_prop_max_nacelle_list = []
P_prop_max_wing_list = []
P_prop_max_fuselage_rear_list = []
P_prop_max_fuselage_underfloor_list = []

# sigma_fcs = sigma_fcs_list_sorted[h]

# Extract input columns from dataframe
index_py = np.vstack(df["index"].apply(ast.literal_eval).to_numpy()) - 1
df["index"] = [tuple(i) for i in index_py]
fcs_loc_unique = df["fcs_loc"].unique()
radius_unique = np.sort(df["radius"].unique())
AR_unique = np.sort(df["AR"].unique())
tdivc_scale_unique = np.sort(df["tdivc_scale"].unique())
N_eng_unique = np.sort(df["N_eng"].unique())
HTR_f_unique = np.sort(df["HTR_f"].unique())
l_fcs_unique = np.sort(df["l_fcs"].unique())
theta_floor_unique = np.sort(df["theta_floor"].unique())
fcs_fuselage_location_unique = np.sort(df["fcs_fuselage_location"].unique())
has_strut_unique = np.sort(df["has_strut"].unique())

# Extract output columns from dataframes

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
P_prop_max_ref = df_ref["P_prop_max"]

# Extract indices of reference aircraft in dataframe columns
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

#%% Data filtering

# Nacelles
mask_nacelle = df["index"].apply(
    lambda x:
    (x[0] == 0) and
    (x[1] == radius_ref_idx) and
    # (x[1] == 19) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_nacelle = df[mask_nacelle]
# print("df_filtered_nacelle =", df_filtered_nacelle)

# Wing
mask_wing = df["index"].apply(
    lambda x:
    (x[0] == 1) and
    (x[1] == radius_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_wing = df[mask_wing]
# print("df_filtered_wing =", df_filtered_wing)

# Fuselage

mask_fuselage_rear = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_fuselage_rear = df[mask_fuselage_rear]
# print("df_filtered_fuselage_rear =", df_filtered_fuselage_rear)

mask_fuselage_underfloor = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    # (x[2] == AR_ref_idx) and
    (x[2] == 8) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[8] == 0) and
    (x[9] == 0)
)
df_filtered_fuselage_underfloor = df[mask_fuselage_underfloor]
# print("df_filtered_fuselage_underfloor =", df_filtered_fuselage_underfloor)

# sys.exit()

#%%

Vol_nacelle = df_filtered_nacelle["Vol_nacelle"].to_numpy()
Vol_wing = df_filtered_wing["Vol_wing"].to_numpy()
Vol_fuselage_rear = df_filtered_fuselage_rear["Vol_fuse"].to_numpy()
Vol_fuselage_underfloor = df_filtered_fuselage_underfloor["Vol_fuse"].to_numpy()

Vol_nacelle_list.append(Vol_nacelle)
Vol_wing_list.append(Vol_wing)
Vol_fuselage_rear_list.append(Vol_fuselage_rear)
Vol_fuselage_underfloor_list.append(Vol_fuselage_underfloor)

PFEI_nacelle = df_filtered_nacelle["PFEI"].to_numpy()
PFEI_wing = df_filtered_wing["PFEI"].to_numpy()
PFEI_fuselage_rear = df_filtered_fuselage_rear["PFEI"].to_numpy()
PFEI_fuselage_underfloor = df_filtered_fuselage_underfloor["PFEI"].to_numpy()

PFEI_nacelle_list.append(PFEI_nacelle)
PFEI_wing_list.append(PFEI_wing)
PFEI_fuselage_rear_list.append(PFEI_fuselage_rear)
PFEI_fuselage_underfloor_list.append(PFEI_fuselage_underfloor)

WMTO_nacelle = df_filtered_nacelle["WMTO"].to_numpy()
WMTO_wing = df_filtered_wing["WMTO"].to_numpy()
WMTO_fuselage_rear = df_filtered_fuselage_rear["WMTO"].to_numpy()
WMTO_fuselage_underfloor = df_filtered_fuselage_underfloor["WMTO"].to_numpy()

WMTO_nacelle_list.append(WMTO_nacelle)
WMTO_wing_list.append(WMTO_wing)
WMTO_fuselage_rear_list.append(WMTO_fuselage_rear)
WMTO_fuselage_underfloor_list.append(WMTO_fuselage_underfloor)

P_prop_max_nacelle = df_filtered_nacelle["P_prop_max"].to_numpy()
P_prop_max_wing = df_filtered_wing["P_prop_max"].to_numpy()
P_prop_max_fuselage_rear = df_filtered_fuselage_rear["P_prop_max"].to_numpy()
P_prop_max_fuselage_underfloor = df_filtered_fuselage_underfloor["P_prop_max"].to_numpy()

P_prop_max_nacelle_list.append(P_prop_max_nacelle)
P_prop_max_wing_list.append(P_prop_max_wing)
P_prop_max_fuselage_rear_list.append(P_prop_max_fuselage_rear)
P_prop_max_fuselage_underfloor_list.append(P_prop_max_fuselage_underfloor)

#%%

import alphashape
from alpha_shape_creator import ensure_xy_array, add_shapely_polygon
from paretoset import paretoset

eta_packing_nacelle_list = [1, 0.55]
eta_packing_wing_list = [1, 0.36]
eta_packing_fuselage_rear_list = [1, 0.69]
eta_packing_fuselage_underfloor_list = [1, 0.88]

nu_req_nacelle = np.concatenate(P_prop_max_nacelle_list) / np.concatenate(Vol_nacelle_list) / 1e6 / eta_packing_nacelle_list[0]
nu_req_wing = np.concatenate(P_prop_max_wing_list) / np.concatenate(Vol_wing_list) / 1e6 / eta_packing_wing_list[0]
nu_req_fuselage_rear = np.concatenate(P_prop_max_fuselage_rear_list) / np.concatenate(Vol_fuselage_rear_list) / 1e6 / eta_packing_fuselage_rear_list[0]
nu_req_fuselage_underfloor = np.concatenate(P_prop_max_fuselage_underfloor_list) / np.concatenate(Vol_fuselage_underfloor_list) / 1e6 / eta_packing_fuselage_underfloor_list[0]

marker_scatter = '.'
alpha_scatter = 0.1
lw_plot = 2

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(nu_req_nacelle, np.concatenate(PFEI_nacelle_list), alpha=alpha_scatter, marker=marker_scatter, edgecolor='none')
ax.scatter(nu_req_wing, np.concatenate(PFEI_wing_list), alpha=alpha_scatter, marker=marker_scatter, edgecolor='none')
ax.scatter(nu_req_fuselage_rear, np.concatenate(PFEI_fuselage_rear_list), alpha=alpha_scatter, marker=marker_scatter, edgecolor='none')
ax.scatter(nu_req_fuselage_underfloor, np.concatenate(PFEI_fuselage_underfloor_list), alpha=alpha_scatter, zorder=-100, marker=marker_scatter, edgecolor='none')

alpha_pareto_list = [0.3, 1]
for packing_idx in [0, 1]:
    
    alpha_pareto = alpha_pareto_list[packing_idx]
    
    nu_req_nacelle = np.concatenate(P_prop_max_nacelle_list) / np.concatenate(Vol_nacelle_list) / 1e6 / eta_packing_nacelle_list[packing_idx]
    nu_req_wing = np.concatenate(P_prop_max_wing_list) / np.concatenate(Vol_wing_list) / 1e6 / eta_packing_wing_list[packing_idx]
    nu_req_fuselage_rear = np.concatenate(P_prop_max_fuselage_rear_list) / np.concatenate(Vol_fuselage_rear_list) / 1e6 / eta_packing_fuselage_rear_list[packing_idx]
    nu_req_fuselage_underfloor = np.concatenate(P_prop_max_fuselage_underfloor_list) / np.concatenate(Vol_fuselage_underfloor_list) / 1e6 / eta_packing_fuselage_underfloor_list[packing_idx]

    # =============================================================================
    points_2d = np.array([
        nu_req_nacelle,
        np.concatenate(PFEI_nacelle_list),
    ]).T
    # pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
    # alpha_shape = alphashape.alphashape(pts, 0.08)
    # add_shapely_polygon(
    #     ax, alpha_shape, facecolor=colors[0], edgecolor='None', alpha=0.5,
    # )
    # cvhull = ConvexHull(points_2d)
    # cvhull_vertices = points_2d[cvhull.vertices]                  # ordered hull vertices
    # polygon = Polygon(cvhull_vertices)
    # x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    # ax.fill(
    #     x_outline_polygon, y_outline_polygon,
    #     hatch='..', facecolor=colors[0], edgecolor='None', alpha=0.3,
    # )
    # ax.plot(x_outline_polygon[:-1], y_outline_polygon[:-1])
    pareto_front = points_2d[paretoset(points_2d, sense=["min", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[0], linewidth=lw_plot, alpha=alpha_pareto)
    pareto_front = points_2d[paretoset(points_2d, sense=["max", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[0], label='Nacelles', linewidth=lw_plot, alpha=alpha_pareto)
    ax.scatter(pareto_front_sorted[:,0][np.argmin(pareto_front_sorted[:,1])], np.min(pareto_front_sorted[:,1]), color=colors[0])
    
    points_2d = np.array([
        nu_req_wing,
        np.concatenate(PFEI_wing_list),
    ]).T
    # pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
    # alpha_shape = alphashape.alphashape(pts, 0.08)
    # add_shapely_polygon(
    #     ax, alpha_shape, facecolor=colors[1], edgecolor='None', alpha=0.5,
    # )
    cvhull = ConvexHull(points_2d)
    cvhull_vertices = points_2d[cvhull.vertices]                  # ordered hull vertices
    polygon = Polygon(cvhull_vertices)
    x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    # ax.fill(
    #     x_outline_polygon, y_outline_polygon,
    #     hatch='..', facecolor=colors[1], edgecolor='None', alpha=0.3,
    # )
    # ax.plot(x_outline_polygon[4:-1], y_outline_polygon[4:-1])
    pareto_front = points_2d[paretoset(points_2d, sense=["min", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[1], linewidth=lw_plot, alpha=alpha_pareto)
    pareto_front = points_2d[paretoset(points_2d, sense=["max", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[1], label='Wing', linewidth=lw_plot, alpha=alpha_pareto)
    ax.scatter(pareto_front_sorted[:,0][np.argmin(pareto_front_sorted[:,1])], np.min(pareto_front_sorted[:,1]), color=colors[1])
    
    points_2d = np.array([
        nu_req_fuselage_rear,
        np.concatenate(PFEI_fuselage_rear_list),
    ]).T
    mask = np.all(np.isfinite(points_2d), axis=1)
    # pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
    # alpha_shape = alphashape.alphashape(pts[mask], 0.01)
    # add_shapely_polygon(
    #     ax, alpha_shape, facecolor=colors[2], edgecolor='None', alpha=0.5,
    # )
    cvhull = ConvexHull(points_2d[mask])
    cvhull_vertices = points_2d[mask][cvhull.vertices]                  # ordered hull vertices
    polygon = Polygon(cvhull_vertices)
    x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    # ax.fill(
    #     x_outline_polygon, y_outline_polygon,
    #     hatch='..', facecolor=colors[2], edgecolor='None', alpha=0.3,
    # )
    pareto_front = points_2d[paretoset(points_2d, sense=["min", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[2], linewidth=lw_plot, alpha=alpha_pareto)
    # ax.scatter(pareto_front_sorted[:,0][np.argmin(pareto_front_sorted[:,1])], np.min(pareto_front_sorted[:,1]), color=colors[2])
    ax.scatter(pareto_front_sorted[:,0][np.argmin(pareto_front_sorted[:,1]) - 1], pareto_front_sorted[:,1][np.argmin(pareto_front_sorted[:,1]) - 1], color=colors[2])
    print(pareto_front_sorted[:,0][np.argmin(pareto_front_sorted[:,1])], np.min(pareto_front_sorted[:,1]))
    print(pareto_front_sorted[:,0], pareto_front_sorted[:,1])
    print()
    pareto_front = points_2d[paretoset(points_2d, sense=["max", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[2], label='Behind cabin', linewidth=lw_plot, alpha=alpha_pareto)
    
    points_2d = np.array([
        nu_req_fuselage_underfloor,
        np.concatenate(PFEI_fuselage_underfloor_list),
    ]).T
    # pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
    # alpha_shape = alphashape.alphashape(pts, 0.1)
    # add_shapely_polygon(
    #     ax, alpha_shape, facecolor=colors[3], edgecolor='None', alpha=0.5,
    # )
    
    cvhull = ConvexHull(points_2d)
    cvhull_vertices = points_2d[cvhull.vertices]                  # ordered hull vertices
    polygon = Polygon(cvhull_vertices)
    x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    # ax.fill(
    #     x_outline_polygon, y_outline_polygon,
    #     hatch='..', facecolor=colors[3], edgecolor='None', alpha=0.3,
    # )
    pareto_front = points_2d[paretoset(points_2d, sense=["min", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    pareto_front_sorted_end = pareto_front_sorted[-1]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[3], linewidth=lw_plot, alpha=alpha_pareto)
    pareto_front = points_2d[paretoset(points_2d, sense=["max", "min"])]
    idx = np.argsort(pareto_front[:, 0])
    pareto_front_sorted = pareto_front[idx]
    pareto_front_sorted_start = pareto_front_sorted[0]
    ax.plot(pareto_front_sorted[:,0], pareto_front_sorted[:,1], color=colors[3], label='Underneath cabin', linewidth=lw_plot, alpha=alpha_pareto)
    ax.plot([pareto_front_sorted_start[0], pareto_front_sorted_end[0]], [pareto_front_sorted_start[1], pareto_front_sorted_end[1]], color=colors[3], linewidth=lw_plot, alpha=alpha_pareto)
    ax.scatter(pareto_front_sorted[:,0][np.argmin(pareto_front_sorted[:,1])], np.min(pareto_front_sorted[:,1]), color=colors[3])
    # =============================================================================

ax.set_xlim(0, 8)

ax.set_xlabel('FCS power density (kW/l)', labelpad=10)
ax.set_ylabel('Relative passenger fuel emission index (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

ax.legend(frameon=False, loc='upper right')
# plt.savefig("PFEI_vs_nu.png", format='png', dpi=600)

plt.show()
