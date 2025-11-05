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

# df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_newest.csv'))
df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_newestttt.csv'))
df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_newest.csv'))

#%%

import ast

# df = sheet_dict['results']

# index = np.array(ast.literal_eval(df["index"].to_numpy(dtype=np.ndarray))) - 1
index_py = np.vstack(df["index"].apply(ast.literal_eval).to_numpy()) - 1
df["index"] = [tuple(i) for i in index_py]
fcs_loc_unique = df["fcs_loc"].unique()
radius_unique = np.sort(df["radius"].unique())
AR_unique = np.sort(df["AR"].unique())
tdivc_scale_unique = np.sort(df["tdivc_scale"].unique())
N_eng_unique = np.sort(df["N_eng"].unique())
HTR_f_unique = np.sort(df["HTR_f"].unique())
Vspec_unique = np.sort(df["Vspec"].unique())
fcs_fuselage_location_unique = np.sort(df["fcs_fuselage_location"].unique())
has_strut_unique = np.sort(df["has_strut"].unique())

# Create 2D coordinate grids
fcs_loc_grid, radius_grid, AR_grid, tdivc_scale_grid, N_eng_grid, HTR_f_grid, Vspec_grid, fcs_fuselage_location_grid, has_strut_grid = np.meshgrid(
    fcs_loc_unique, radius_unique, AR_unique, tdivc_scale_unique, N_eng_unique, HTR_f_unique, Vspec_unique, fcs_fuselage_location_unique, has_strut_unique, indexing="ij"
)

# Create shape tuple for reshaping
shape = (
    len(fcs_loc_unique),
    len(radius_unique),
    len(AR_unique),
    len(tdivc_scale_unique),
    len(N_eng_unique),
    len(HTR_f_unique),
    len(Vspec_unique),
    len(fcs_fuselage_location_unique),
    len(has_strut_unique),
    
)

# Create index mapping for each dimension
fcs_loc_idx = {v: i for i, v in enumerate(fcs_loc_unique)}
radius_idx = {v: i for i, v in enumerate(radius_unique)}
AR_idx = {v: i for i, v in enumerate(AR_unique)}
tdivc_scale_idx = {v: i for i, v in enumerate(tdivc_scale_unique)}
N_eng_idx = {v: i for i, v in enumerate(N_eng_unique)}
HTR_f_idx = {v: i for i, v in enumerate(HTR_f_unique)}
Vspec_idx = {v: i for i, v in enumerate(Vspec_unique)}
fcs_fuselage_location_idx = {v: i for i, v in enumerate(fcs_fuselage_location_unique)}
has_strut_idx = {v: i for i, v in enumerate(has_strut_unique)}

# Initialize empty arrays
dtype = np.float32
CDS_grid = np.full(shape, np.nan, dtype=dtype)
WMTO_grid = np.full(shape, np.nan, dtype=dtype)
Vol_wing_grid = np.full(shape, np.nan, dtype=dtype)
PFEI_grid = np.full(shape, np.nan, dtype=dtype)
seats_abreast_grid = np.full(shape, np.nan, dtype=dtype)
Vol_nacelle_grid = np.full(shape, np.nan, dtype=dtype)
y_centroid_wing_grid = np.full(shape, np.nan, dtype=dtype)
L_fuse_grid = np.full(shape, np.nan, dtype=dtype)
y_centroid_nacelles_grid = np.full(shape, np.nan, dtype=dtype)
span_grid = np.full(shape, np.nan, dtype=dtype)

# Populate arrays
for _, row in df.iterrows():
    h = fcs_loc_idx[row["fcs_loc"]]
    i = radius_idx[row["radius"]]
    j = AR_idx[row["AR"]]
    k = tdivc_scale_idx[row["tdivc_scale"]]
    l = N_eng_idx[row["N_eng"]]
    m = HTR_f_idx[row["HTR_f"]]
    n = Vspec_idx[row["Vspec"]]
    o = fcs_fuselage_location_idx[row["fcs_fuselage_location"]]
    p = has_strut_idx[row["has_strut"]]
    
    CDS_grid[h, i, j, k, l, m, n, o, p] = row["CDS"]
    WMTO_grid[h, i, j, k, l, m, n, o, p] = row["WMTO"]
    Vol_wing_grid[h, i, j, k, l, m, n, o, p] = row["Vol_wing"]
    PFEI_grid[h, i, j, k, l, m, n, o, p] = row["PFEI"]
    seats_abreast_grid[h, i, j, k, l, m, n, o, p] = row["seats_abreast"]
    Vol_nacelle_grid[h, i, j, k, l, m, n, o, p] = row["Vol_nacelle"]
    y_centroid_wing_grid[h, i, j, k, l, m, n, o, p] = row["y_centroid_wing"]
    L_fuse_grid[h, i, j, k, l, m, n, o, p] = row["L_fuse"]
    y_centroid_nacelles_grid[h, i, j, k, l, m, n, o, p] = row["y_centroid_nacelles"]
    span_grid[h, i, j, k, l, m, n, o, p] = row["span"]
    
sys.exit()

y_centroid_wing_norm_grid = y_centroid_wing_grid / (span_grid / 2)
y_centroid_nacelles_norm_grid = y_centroid_nacelles_grid / (span_grid / 2)

#%%

fcs_loc_ref = df_ref["fcs_loc"]
radius_ref = df_ref["radius"]
AR_ref = df_ref["AR"]
tdivc_scale_ref = df_ref["tdivc_scale"]
N_eng_ref = df_ref["N_eng"]
HTR_f_ref = df_ref["HTR_f"]
Vspec_ref = df_ref["Vspec"]
fcs_fuselage_location_ref = df_ref["fcs_fuselage_location"]
has_strut_ref = df_ref["has_strut"]

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

#%%

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









