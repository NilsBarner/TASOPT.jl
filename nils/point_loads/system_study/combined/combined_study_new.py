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
# # df = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'volume_results_narrowbody_nacelle_newest.csv'))
# df = pd.read_csv(os.path.join(os.getcwd(), 'nils', 'point_loads', 'system_study', 'volume', 'volume_results_narrowbody_newesttttt.csv'))
df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_newest_1500.csv'))
df = pd.read_csv(os.path.join(os.getcwd(), 'combined_results_narrowbody_1500.csv'))
# df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_newest_3000.csv'))
# df = pd.read_csv(os.path.join(os.getcwd(), 'combined_results_narrowbody_3000.csv'))

# =============================================================================
#%%
import os
import glob
import re
import pandas as pd

folder = os.getcwd()  # or any other folder path

# Find all CSV files that start with 'combined_' and end with '_<number>.csv'
csv_files = glob.glob(os.path.join(folder, 'combined_*_[0-9]*.csv'))

# Sort numerically by the trailing number (so ..._2.csv comes before ..._10.csv)
def extract_number(filename):
    match = re.search(r'_(\d+)\.csv$', filename)
    return int(match.group(1)) if match else -1

csv_files = sorted(csv_files, key=extract_number)
sigma_fcs_list_sorted = [extract_number(csv_file) for csv_file in csv_files]

# Read all into a list of DataFrames
dfs = [pd.read_csv(f) for f in csv_files]

print(f"Loaded {len(dfs)} matching CSV files.")

#%%
# =============================================================================

nu_fcs_range = np.linspace(0.2, 5.5, 1000) * 1e6
# nu_fcs_range = [5.5e6]
# nu_fcs_range = [1e6]

X, Y = np.meshgrid(nu_fcs_range, sigma_fcs_list_sorted)

PFEI_min_nacelle_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
PFEI_min_wing_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
PFEI_min_fuselage_rear_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
PFEI_min_fuselage_underfloor_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan

WMTO_min_nacelle_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
WMTO_min_wing_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
WMTO_min_fuselage_rear_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
WMTO_min_fuselage_underfloor_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan

CDS_min_nacelle_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
CDS_min_wing_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
CDS_min_fuselage_rear_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan
CDS_min_fuselage_underfloor_array = np.ones((len(dfs), len(nu_fcs_range))) * np.nan

for h, df in enumerate(dfs):
    
    sigma_fcs = sigma_fcs_list_sorted[h]

    # Extract input columns from dataframe
    index_py = np.vstack(df["index"].apply(ast.literal_eval).to_numpy()) - 1
    df["index"] = [tuple(i) for i in index_py]
    sigma_fcs_unique = df["sigma_fcs"].unique()
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
    
    for i, nu_fcs in enumerate(nu_fcs_range):
    
        # Vol_fcs_nacelle_req = df_filtered_nacelle['P_prop_max'] / nu_fcs
        # Vol_fcs_wing_req = df_filtered_wing['P_prop_max'] / nu_fcs
        # Vol_fcs_fuselage_rear_req = df_filtered_fuselage_rear['P_prop_max'] / nu_fcs
        # Vol_fcs_fuselage_underfloor_req = df_filtered_fuselage_underfloor['P_prop_max'] / nu_fcs
    
        # # Filter out all designs where Vol_fcs_req < Vol_fcs_av
        
        # mask_nacelle_vol = df_filtered_nacelle["Vol_nacelle"].apply(lambda x: x > Vol_fcs_nacelle_req)
        # df_filtered_nacelle_vol = df_filtered_nacelle[mask_nacelle_vol]
        # mask_wing_vol = df_filtered_wing["Vol_wing"].apply(lambda x: x > Vol_fcs_wing_req)
        # df_filtered_wing_vol = df_filtered_wing[mask_wing_vol]
        
        # sys.exit('Stop.')
        # mask_fuselage_rear_vol = df_filtered_fuselage_rear["Vol_fuse"].apply(lambda x: x > Vol_fcs_fuselage_rear_req)
        # df_filtered_fuselage_rear_vol = df_filtered_fuselage_rear[mask_fuselage_rear_vol]
        # mask_fuselage_underfloor_vol = df_filtered_fuselage_underfloor["Vol_fuse"].apply(lambda x: x > Vol_fcs_fuselage_underfloor_req)
        # df_filtered_fuselage_underfloor_vol = df_filtered_fuselage_underfloor[mask_fuselage_underfloor_vol]
        
        Vol_fcs_nacelle_req = df_filtered_nacelle['P_prop_max'] / nu_fcs
        Vol_fcs_wing_req = df_filtered_wing['P_prop_max'] / nu_fcs
        Vol_fcs_fuselage_rear_req = df_filtered_fuselage_rear['P_prop_max'] / nu_fcs
        Vol_fcs_fuselage_underfloor_req = df_filtered_fuselage_underfloor['P_prop_max'] / nu_fcs
        
        # correct masks (returns a 1-D boolean Series)
        mask_nacelle_vol = df_filtered_nacelle["Vol_nacelle"] > Vol_fcs_nacelle_req
        mask_wing_vol = df_filtered_wing["Vol_wing"] > Vol_fcs_wing_req
        mask_fuselage_rear_vol = df_filtered_fuselage_rear["Vol_fuse"] > Vol_fcs_fuselage_rear_req
        mask_fuselage_underfloor_vol = df_filtered_fuselage_underfloor["Vol_fuse"] > Vol_fcs_fuselage_underfloor_req
        
        df_filtered_nacelle_vol = df_filtered_nacelle[mask_nacelle_vol]
        df_filtered_wing_vol = df_filtered_wing[mask_wing_vol]
        # print('np.shape(df_filtered_wing_vol), np.shape(df_filtered_wing) =', np.shape(df_filtered_wing_vol), np.shape(df_filtered_wing))
        df_filtered_fuselage_rear_vol = df_filtered_fuselage_rear[mask_fuselage_rear_vol]
        df_filtered_fuselage_underfloor_vol = df_filtered_fuselage_underfloor[mask_fuselage_underfloor_vol]
    
        PFEI_min_nacelle = min(df_filtered_nacelle_vol["PFEI"]) if len(df_filtered_nacelle_vol["PFEI"]) > 0 else np.nan
        PFEI_min_wing = min(df_filtered_wing_vol["PFEI"]) if len(df_filtered_wing_vol["PFEI"]) > 0 else np.nan
        PFEI_min_fuselage_rear = min(df_filtered_fuselage_rear_vol["PFEI"]) if len(df_filtered_fuselage_rear_vol["PFEI"]) > 0 else np.nan
        PFEI_min_fuselage_underfloor = min(df_filtered_fuselage_underfloor_vol["PFEI"]) if len(df_filtered_fuselage_underfloor_vol["PFEI"]) > 0 else np.nan
        
        WMTO_min_nacelle = min(df_filtered_nacelle_vol["WMTO"]) if len(df_filtered_nacelle_vol["WMTO"]) > 0 else np.nan
        WMTO_min_wing = min(df_filtered_wing_vol["WMTO"]) if len(df_filtered_wing_vol["WMTO"]) > 0 else np.nan
        WMTO_min_fuselage_rear = min(df_filtered_fuselage_rear_vol["WMTO"]) if len(df_filtered_fuselage_rear_vol["WMTO"]) > 0 else np.nan
        WMTO_min_fuselage_underfloor = min(df_filtered_fuselage_underfloor_vol["WMTO"]) if len(df_filtered_fuselage_underfloor_vol["WMTO"]) > 0 else np.nan
        
        CDS_min_nacelle = min(df_filtered_nacelle_vol["CDS"]) if len(df_filtered_nacelle_vol["CDS"]) > 0 else np.nan
        CDS_min_wing = min(df_filtered_wing_vol["CDS"]) if len(df_filtered_wing_vol["CDS"]) > 0 else np.nan
        CDS_min_fuselage_rear = min(df_filtered_fuselage_rear_vol["CDS"]) if len(df_filtered_fuselage_rear_vol["CDS"]) > 0 else np.nan
        CDS_min_fuselage_underfloor = min(df_filtered_fuselage_underfloor_vol["CDS"]) if len(df_filtered_fuselage_underfloor_vol["CDS"]) > 0 else np.nan
        
        PFEI_min_nacelle_array[h, i] = PFEI_min_nacelle
        PFEI_min_wing_array[h, i] = PFEI_min_wing
        PFEI_min_fuselage_rear_array[h, i] = PFEI_min_fuselage_rear
        PFEI_min_fuselage_underfloor_array[h, i] = PFEI_min_fuselage_underfloor
        
        WMTO_min_nacelle_array[h, i] = WMTO_min_nacelle
        WMTO_min_wing_array[h, i] = WMTO_min_wing
        WMTO_min_fuselage_rear_array[h, i] = WMTO_min_fuselage_rear
        WMTO_min_fuselage_underfloor_array[h, i] = WMTO_min_fuselage_underfloor
        
        CDS_min_nacelle_array[h, i] = CDS_min_nacelle
        CDS_min_wing_array[h, i] = CDS_min_wing
        CDS_min_fuselage_rear_array[h, i] = CDS_min_fuselage_rear
        CDS_min_fuselage_underfloor_array[h, i] = CDS_min_fuselage_underfloor
    
    ###

    # ax.plot(nu_fcs_range, PFEI_min_nacelle_list, color=colors[h])
    # ax.plot(nu_fcs_range, PFEI_min_wing_list, color=colors[h])
    # ax.plot(nu_fcs_range, PFEI_min_fuselage_rear_list, color=colors[h])
    # ax.plot(nu_fcs_range, PFEI_min_fuselage_underfloor_list, color=colors[h])

#%%

stacked = np.stack([
        PFEI_min_nacelle_array,
        PFEI_min_wing_array,
        PFEI_min_fuselage_rear_array,
        PFEI_min_fuselage_underfloor_array,
    ], axis=0
)

# stacked = np.stack([
#         WMTO_min_nacelle_array,
#         WMTO_min_wing_array,
#         WMTO_min_fuselage_rear_array,
#         WMTO_min_fuselage_underfloor_array,
#     ], axis=0
# )

# stacked = np.stack([
#         CDS_min_nacelle_array,
#         CDS_min_wing_array,
#         CDS_min_fuselage_rear_array,
#         CDS_min_fuselage_underfloor_array,
#     ], axis=0
# )

# Compute elementwise minimum ignoring NaNs
min_vals = np.nanmin(stacked, axis=0)

# Get indices (0–3) of the array each minimum came from
min_idx = np.nanargmin(stacked, axis=0)

# # =============================================================================
# # Compute elementwise minimum ignoring NaNs
# min_vals = np.nanmax(stacked, axis=0)

# # Get indices (0–3) of the array each minimum came from
# min_idx = np.nanargmax(stacked, axis=0)
# # =============================================================================

from matplotlib.colors import ListedColormap, BoundaryNorm

# Create a discrete colormap
cmap = ListedColormap(colors[:4])

# Define boundaries between color bins
bounds = np.arange(-0.5, 4.5, 1)  # [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
norm = BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(8, 6))

# pcm = ax.pcolormesh(X, Y, PFEI_min_nacelle_array, cmap='viridis')
ax.imshow(
    min_idx, extent=[X.min() / 1e6, X.max() / 1e6, Y.min() / 1e3, Y.max() / 1e3],
    cmap=cmap, norm=norm, origin='lower', aspect='auto',
)
ax.set_xlabel('FCS power density (kW/l)', labelpad=10)
ax.set_ylabel('FCS specific power (kW/kg)', labelpad=10)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)
ax.grid(True)

legend_elements = [
    Patch(facecolor=colors[0], edgecolor='None', alpha=1, label='Nacelle'),
    Patch(facecolor=colors[1], edgecolor='None', alpha=1, label='Wing'),
    Patch(facecolor=colors[2], edgecolor='None', alpha=1, label='Behind cabin'),
    Patch(facecolor=colors[3], edgecolor='None', alpha=1, label='Underneath cabin'),
]
leg = ax.legend(handles=legend_elements, loc='upper right', frameon=True)

plt.show()

# fig, ax = plt.subplots()
# plt.show()


