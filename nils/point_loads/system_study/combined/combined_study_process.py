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

# Import data from .csv files
# df = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'volume_results_narrowbody_nacelle_newest.csv'))
df = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'volume_results_narrowbody_newesttttt.csv'))
df_ref = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'volume_results_narrowbody_ref_newest.csv'))

# Extract input columns from dataframe
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

#%% Data filtering

# Nacelles
mask_nacelle = df["index"].apply(
    lambda x:
    (x[0] == 0) and
    (x[1] == radius_ref_idx) and
    (x[2] == AR_ref_idx) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[6] == l_fcs_ref_idx) and
    (x[7] == theta_floor_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_nacelle = df[mask_nacelle]

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

# Fuselage

# mask_fuselage = df["index"].apply(
#     lambda x:
#     (x[0] == 2) and
#     (x[2] == AR_ref_idx) and
#     (x[3] == tdivc_scale_ref_idx) and
#     (x[4] == N_eng_ref_idx) and
#     (x[5] == HTR_f_ref_idx) and
#     (x[9] == 0)
# )
# df_filtered_fuselage = df[mask_fuselage]

mask_fuselage_rear = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    (x[2] == AR_ref_idx) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[8] == 1) and
    (x[9] == 0)
)
df_filtered_fuselage_rear = df[mask_fuselage_rear]

mask_fuselage_underfloor = df["index"].apply(
    lambda x:
    (x[0] == 2) and
    (x[2] == AR_ref_idx) and
    (x[3] == tdivc_scale_ref_idx) and
    (x[4] == N_eng_ref_idx) and
    (x[5] == HTR_f_ref_idx) and
    (x[8] == 0) and
    (x[9] == 0)
)
df_filtered_fuselage_underfloor = df[mask_fuselage_underfloor]

#%%

def extract_loci_from_contour(cs):
    """
    Robustly extract contour segments from a QuadContourSet `cs`.
    Returns a list of arrays, each of shape (M,2) with columns [x, y].
    Works whether cs has .collections or only .allsegs.
    """
    loci = []
    # Preferred: use .allsegs (list-of-levels -> list-of-segs)
    if hasattr(cs, "allsegs") and cs.allsegs:
        # cs.allsegs is a list (one element per level). We'll flatten all segments.
        for level_segs in cs.allsegs:
            for seg in level_segs:
                if seg.size and seg.shape[0] >= 2:
                    loci.append(np.asarray(seg))
    else:
        # Fallback: look for collections then get_paths
        if hasattr(cs, "collections"):
            for coll in cs.collections:
                for path in coll.get_paths():
                    coords = path.vertices
                    if coords.shape[0] >= 2:
                        loci.append(np.asarray(coords))
    return loci


def resample_curve_by_arclength(xy, n_points, ax):
    
    # Transform to display (pixel) coordinates
    disp = ax.transData.transform(xy)
    d = np.hypot(np.diff(disp[:,0]), np.diff(disp[:,1]))
    s = np.concatenate(([0.0], d.cumsum()))
    
    t = np.linspace(0, s[-1], n_points)
    xi = np.interp(t, s, disp[:,0])
    yi = np.interp(t, s, disp[:,1])
    disp_i = np.column_stack((xi, yi))
    
    # Back to data coords
    pts = ax.transData.inverted().transform(disp_i)
    return pts


def resample_all_loci(loci, n_points, ax):
    """Apply the endpoint-normalised resampling to a list of loci."""
    return [resample_curve_by_arclength(locus, n_points, ax) for locus in loci]

#%%

P_fcs = 20e6
# nu_fcs_range = np.linspace(0.2, 1.0, 10) * 1e6
# nu_fcs_range = np.linspace(0.2, 5.5, 10) * 1e6
nu_fcs_range = np.linspace(0.2, 2, 10) * 1e6

# fig, ax = plt.subplots()
# ax.scatter(df_filtered_wing['Vol_wing'], df_filtered_wing['y_centroid_wing'])
# plt.show()
# sys.exit()

df = df_filtered_wing
input_str_1 = 'AR'
input_str_2 = 'tdivc_scale'
output_str = 'Vol_wing'
target = 25

# df = df_filtered_nacelle
# input_str_1 = 'N_eng'
# input_str_2 = 'HTR_f'
# output_str = 'Vol_nacelle'
# target = 25

# df = df_filtered_fuselage_rear
# input_str_1 = 'radius'
# input_str_2 = 'l_fcs'
# output_str = 'Vol_fuse'
# target = 25

# df = df_filtered_fuselage_underfloor
# input_str_1 = 'radius'
# input_str_2 = 'theta_floor'
# output_str = 'Vol_fuse'
# target = 75

"""
location_list = ["wing", "nacelle", "underfloor", "rear"]
input_strs_list = [
    ['AR', 'tdivc_scale'], ['N_eng', 'HTR_f'], ['radius', 'l_fcs'], ['radius', 'theta_floor']
]
input_strs_dict = dict(zip(location_list, input_strs_list))
for location in location_list:
    input_strs = input_strs_dict[location]
    
    
"""

#%%
'''
# Group and unstack: rows = AR, columns = t_c
pivot = df.groupby([input_str_1, input_str_2])[output_str].mean().unstack()

# Correct orientation
AR_vals = pivot.index.values       # x-axis
tc_vals = pivot.columns.values     # y-axis
Z = pivot.values.T                 # transpose to match meshgrid orientation

# Build mesh so x=AR, y=t_c
ARg, tcg = np.meshgrid(AR_vals, tc_vals)

# Plot contour
fig, ax = plt.subplots()
cs = ax.contour(ARg, tcg, Z, levels=[target])

loci = extract_loci_from_contour(cs)
samples = resample_all_loci(loci, n_points=10, ax=ax)
for sample in samples:
    ax.scatter(sample[:,0], sample[:,1])

ax.clabel(cs, inline=True)
ax.set_xlabel('Aspect ratio AR')
ax.set_ylabel('Thickness-to-chord t/c')
plt.show()
'''
#%%

# Sample points and values
pts = df[[input_str_1, input_str_2]].values
vals = df[output_str].values

# Create uniform grid in x=AR, y=t_c
ar_lin = np.linspace(df[input_str_1].min(), df[input_str_1].max(), 200)
tc_lin = np.linspace(df[input_str_2].min(), df[input_str_2].max(), 200)
ARg, tcg = np.meshgrid(ar_lin, tc_lin)

# Interpolate
Zgrid = griddata(pts, vals, (ARg, tcg), method='cubic')

fig, ax = plt.subplots()

np.array

for nu_fcs in nu_fcs_range:
    
    V_fcs = P_fcs / nu_fcs
    cs = ax.contour(ARg, tcg, Zgrid, levels=[V_fcs], colors='r')
    
    loci = extract_loci_from_contour(cs)
    samples = resample_all_loci(loci, n_points=10, ax=ax)
    print('samples =', samples)
    for sample in samples:
        ax.scatter(sample[:,0], sample[:,1])
    
    ax.scatter(df[input_str_1], df[input_str_2], s=5, alpha=0.2)
    
ax.set_xlabel('Aspect ratio AR')
ax.set_ylabel('Thickness-to-chord t/c')
plt.show()

# # Extract contour coordinates (AR, t_c)
# loci = [path.vertices for col in cs.collections for path in col.get_paths()]

#%%

# data = np.array([(1, 2, 3), (4, 5, 6), (7, 8, 9)],
#                 dtype=[("a", "i4"), ("b", "i4"), ("c", "i4")])
# wing_ = pd.DataFrame(data, columns=['c', 'a'])
