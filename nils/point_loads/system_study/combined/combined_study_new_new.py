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
# df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_newest_1500.csv'))
# df = pd.read_csv(os.path.join(os.getcwd(), 'combined_results_narrowbody_1500.csv'))
# df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_newest_3000.csv'))
# df = pd.read_csv(os.path.join(os.getcwd(), 'combined_results_narrowbody_3000.csv'))

df_ref = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_ref_101125_1500.csv'))

# =============================================================================
#%%
import os
import glob
import re
import pandas as pd

folder = os.getcwd()  # or any other folder path

# Find all CSV files that start with 'combined_' and end with '_<number>.csv'
csv_files = glob.glob(os.path.join(folder, 'combined_*narrowbody_[0-9]*.csv'))

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

nu_fcs_range = np.linspace(0.2, 5.5, 10) * 1e6
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
    
    sigma_fcs_nacelle = np.ones_like(df_filtered_nacelle["Vol_nacelle"]) * sigma_fcs
    sigma_fcs_wing = np.ones_like(df_filtered_wing["Vol_wing"]) * sigma_fcs
    sigma_fcs_fuselage_rear = np.ones_like(df_filtered_fuselage_rear["Vol_fuse"]) * sigma_fcs
    sigma_fcs_fuselage_underfloor = np.ones_like(df_filtered_fuselage_underfloor["Vol_fuse"]) * sigma_fcs
    
    sigma_fcs_nacelle_list.append(sigma_fcs_nacelle)
    sigma_fcs_wing_list.append(sigma_fcs_wing)
    sigma_fcs_fuselage_rear_list.append(sigma_fcs_fuselage_rear)
    sigma_fcs_fuselage_underfloor_list.append(sigma_fcs_fuselage_underfloor)
    
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
    
    for i, nu_fcs in enumerate(nu_fcs_range):
    
        nu_fcs_nacelle_req = df_filtered_nacelle['P_prop_max'].to_numpy() / df_filtered_nacelle["Vol_nacelle"].to_numpy()
        nu_fcs_wing_req = df_filtered_wing['P_prop_max'].to_numpy() / df_filtered_wing["Vol_wing"].to_numpy()
        nu_fcs_fuselage_rear_req = df_filtered_fuselage_rear['P_prop_max'].to_numpy() / df_filtered_fuselage_rear["Vol_fuse"].to_numpy()
        nu_fcs_fuselage_underfloor_req = df_filtered_fuselage_underfloor['P_prop_max'].to_numpy() / df_filtered_fuselage_underfloor["Vol_fuse"].to_numpy()
        
        mask_nacelle_vol = nu_fcs > nu_fcs_nacelle_req
        mask_wing_vol = nu_fcs > nu_fcs_wing_req
        mask_fuselage_rear_vol = nu_fcs > nu_fcs_fuselage_rear_req
        mask_fuselage_underfloor_vol = nu_fcs > nu_fcs_fuselage_underfloor_req
    
        # # Filter out all designs where Vol_fcs_req < Vol_fcs_av
        
        # mask_nacelle_vol = df_filtered_nacelle["Vol_nacelle"] > Vol_fcs_nacelle_req
        # mask_wing_vol = df_filtered_wing["Vol_wing"] > Vol_fcs_wing_req
        # mask_fuselage_rear_vol = df_filtered_fuselage_rear["Vol_fuse"] > Vol_fcs_fuselage_rear_req
        # mask_fuselage_underfloor_vol = df_filtered_fuselage_underfloor["Vol_fuse"] > Vol_fcs_fuselage_underfloor_req
        
        df_filtered_nacelle_vol = df_filtered_nacelle[mask_nacelle_vol]
        df_filtered_wing_vol = df_filtered_wing[mask_wing_vol]
        # print('np.shape(df_filtered_wing_vol), np.shape(df_filtered_wing) =', np.shape(df_filtered_wing_vol), np.shape(df_filtered_wing))
        df_filtered_fuselage_rear_vol = df_filtered_fuselage_rear[mask_fuselage_rear_vol]
        df_filtered_fuselage_underfloor_vol = df_filtered_fuselage_underfloor[mask_fuselage_underfloor_vol]
        
        # # =============================================================================
        # df_filtered_nacelle_vol = df_filtered_nacelle
        # df_filtered_wing_vol = df_filtered_wing
        # df_filtered_fuselage_rear_vol = df_filtered_fuselage_rear
        # df_filtered_fuselage_underfloor_vol = df_filtered_fuselage_underfloor
        # # =============================================================================
    
        PFEI_min_nacelle = np.nanmin(df_filtered_nacelle_vol["PFEI"]) if len(df_filtered_nacelle_vol["PFEI"]) > 0 else np.nan
        PFEI_min_wing = np.nanmin(df_filtered_wing_vol["PFEI"]) if len(df_filtered_wing_vol["PFEI"]) > 0 else np.nan
        PFEI_min_fuselage_rear = np.nanmin(df_filtered_fuselage_rear_vol["PFEI"]) if len(df_filtered_fuselage_rear_vol["PFEI"]) > 0 else np.nan
        PFEI_min_fuselage_underfloor = np.nanmin(df_filtered_fuselage_underfloor_vol["PFEI"]) if len(df_filtered_fuselage_underfloor_vol["PFEI"]) > 0 else np.nan
        
        # print(PFEI_min_nacelle, PFEI_min_wing, PFEI_min_fuselage_rear, PFEI_min_fuselage_underfloor)
        # print()
        
        WMTO_min_nacelle = np.nanmin(df_filtered_nacelle_vol["WMTO"]) if len(df_filtered_nacelle_vol["WMTO"]) > 0 else np.nan
        WMTO_min_wing = np.nanmin(df_filtered_wing_vol["WMTO"]) if len(df_filtered_wing_vol["WMTO"]) > 0 else np.nan
        WMTO_min_fuselage_rear = np.nanmin(df_filtered_fuselage_rear_vol["WMTO"]) if len(df_filtered_fuselage_rear_vol["WMTO"]) > 0 else np.nan
        WMTO_min_fuselage_underfloor = np.nanmin(df_filtered_fuselage_underfloor_vol["WMTO"]) if len(df_filtered_fuselage_underfloor_vol["WMTO"]) > 0 else np.nan
        
        CDS_min_nacelle = np.nanmin(df_filtered_nacelle_vol["CDS"]) if len(df_filtered_nacelle_vol["CDS"]) > 0 else np.nan
        CDS_min_wing = np.nanmin(df_filtered_wing_vol["CDS"]) if len(df_filtered_wing_vol["CDS"]) > 0 else np.nan
        CDS_min_fuselage_rear = np.nanmin(df_filtered_fuselage_rear_vol["CDS"]) if len(df_filtered_fuselage_rear_vol["CDS"]) > 0 else np.nan
        CDS_min_fuselage_underfloor = np.nanmin(df_filtered_fuselage_underfloor_vol["CDS"]) if len(df_filtered_fuselage_underfloor_vol["CDS"]) > 0 else np.nan
        
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

# sys.exit()
# '''
#%%

fig, ax = plt.subplots()

# ax.imshow(
#     PFEI_min_wing_array, extent=[X.min() / 1e6, X.max() / 1e6, Y.min() / 1e3, Y.max() / 1e3],
#     # cmap=cmap, norm=norm, origin='lower', aspect='auto',
# )

# Vol_nacelle_list = []
# Vol_wing_list = []
# Vol_fuselage_rear_list = []
# Vol_fuselage_underfloor_list = []

# PFEI_nacelle_list = []
# PFEI_wing_list = []
# PFEI_fuselage_rear_list = []
# PFEI_fuselage_underfloor_list = []

# ax.scatter(np.concatenate(Vol_nacelle_list) / 0.55, np.concatenate(PFEI_nacelle_list), c=np.concatenate(sigma_fcs_nacelle_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_wing_list) / 0.36, np.concatenate(PFEI_wing_list), c=np.concatenate(sigma_fcs_wing_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_fuselage_rear_list) / 0.69, np.concatenate(PFEI_fuselage_rear_list), c=np.concatenate(sigma_fcs_fuselage_rear_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_fuselage_underfloor_list) / 0.88, np.concatenate(PFEI_fuselage_underfloor_list), alpha=0.1, c=np.concatenate(sigma_fcs_fuselage_underfloor_list), zorder=-100)

ax.scatter(np.concatenate(P_prop_max_nacelle_list) / np.concatenate(Vol_nacelle_list) / 0.55 / 1e6, np.concatenate(PFEI_nacelle_list), alpha=0.1)
ax.scatter(np.concatenate(P_prop_max_wing_list) / np.concatenate(Vol_wing_list) / 0.36 / 1e6, np.concatenate(PFEI_wing_list), alpha=0.1)
ax.scatter(np.concatenate(P_prop_max_fuselage_rear_list) / np.concatenate(Vol_fuselage_rear_list) / 0.69 / 1e6, np.concatenate(PFEI_fuselage_rear_list), alpha=0.1)
ax.scatter(np.concatenate(P_prop_max_fuselage_underfloor_list) / np.concatenate(Vol_fuselage_underfloor_list) / 0.88 / 1e6, np.concatenate(PFEI_fuselage_underfloor_list) * 1.1, alpha=0.1, zorder=-100)

# ax.scatter(np.concatenate(P_prop_max_nacelle_list) / np.concatenate(Vol_nacelle_list) / 0.55 / 1e6, np.concatenate(PFEI_nacelle_list), c=np.concatenate(sigma_fcs_nacelle_list), alpha=0.1)
# ax.scatter(np.concatenate(P_prop_max_wing_list) / np.concatenate(Vol_wing_list) / 0.36 / 1e6, np.concatenate(PFEI_wing_list), c=np.concatenate(sigma_fcs_wing_list), alpha=0.1)
# ax.scatter(np.concatenate(P_prop_max_fuselage_rear_list) / np.concatenate(Vol_fuselage_rear_list) / 0.69 / 1e6, np.concatenate(PFEI_fuselage_rear_list), c=np.concatenate(sigma_fcs_fuselage_rear_list), alpha=0.1)
# ax.scatter(np.concatenate(P_prop_max_fuselage_underfloor_list) / np.concatenate(Vol_fuselage_underfloor_list) / 0.88 / 1e6, np.concatenate(PFEI_fuselage_underfloor_list), c=np.concatenate(sigma_fcs_fuselage_underfloor_list), alpha=0.1, zorder=-100)

# ax.scatter(np.concatenate(Vol_nacelle_list) / 0.55, np.concatenate(P_prop_max_nacelle_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_wing_list) / 0.36, np.concatenate(P_prop_max_wing_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_fuselage_rear_list) / 0.69, np.concatenate(P_prop_max_fuselage_rear_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_fuselage_underfloor_list) / 0.88, np.concatenate(P_prop_max_fuselage_underfloor_list), alpha=0.1, zorder=-100)

# ax.scatter(np.concatenate(Vol_nacelle_list) / 0.55, np.concatenate(WMTO_nacelle_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_wing_list) / 0.36, np.concatenate(WMTO_wing_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_fuselage_rear_list) / 0.69, np.concatenate(WMTO_fuselage_rear_list), alpha=0.1)
# ax.scatter(np.concatenate(Vol_fuselage_underfloor_list) / 0.88, np.concatenate(WMTO_fuselage_underfloor_list), alpha=0.1, zorder=-100)

# for (Vol_nacelle, PFEI_nacelle) in zip(Vol_nacelle_list, PFEI_nacelle_list):
#     ax.scatter(Vol_nacelle / 0.55, PFEI_nacelle)
    
# for (Vol_wing, PFEI_wing) in zip(Vol_wing_list, PFEI_wing_list):
#     ax.scatter(Vol_wing / 0.36, PFEI_wing)
    
# for (Vol_fuselage_rear, PFEI_fuselage_rear) in zip(Vol_fuselage_rear_list, PFEI_fuselage_rear_list):
#     ax.scatter(Vol_fuselage_rear / 0.55, PFEI_fuselage_rear)
    
# for (Vol_fuselage_underfloor, PFEI_fuselage_underfloor) in zip(Vol_fuselage_underfloor_list, PFEI_fuselage_underfloor_list):
#     ax.scatter(Vol_fuselage_underfloor / 0.55, PFEI_fuselage_underfloor)

ax.set_xlim(0, 10)

plt.show()

sys.exit()

#%%
# '''
stacked = np.stack([
        PFEI_min_nacelle_array,
        PFEI_min_wing_array,
        PFEI_min_fuselage_rear_array,
        PFEI_min_fuselage_underfloor_array * 1.15,
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

# Mask positions where all arrays are NaN
all_nan_mask = np.all(np.isnan(stacked), axis=0)

# For argmin, temporarily replace NaNs with +inf so argmin skips them safely
safe_stacked = np.where(np.isnan(stacked), np.inf, stacked)
min_idx = np.argmin(safe_stacked, axis=0)

# # Compute elementwise minimum ignoring NaNs
# min_vals = np.nanmin(stacked, axis=0)

# # Get indices (0–3) of the array each minimum came from
# min_idx = np.nanargmin(stacked, axis=0)

# Assign a sentinel (e.g. -1) where all values were NaN
min_idx[all_nan_mask] = -1

# # =============================================================================
# # Compute elementwise minimum ignoring NaNs
# min_vals = np.nanmax(stacked, axis=0)

# # Get indices (0–3) of the array each minimum came from
# min_idx = np.nanargmax(stacked, axis=0)
# # =============================================================================

from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.ndimage import zoom

# --- upsample (smooth) the label field to a finer grid ---
upsample_factor = 10              # larger -> smoother boundaries
# zoom with spline interpolation (order=3) -> smooth continuous field
smooth_field = zoom(min_idx.astype(float), zoom=upsample_factor, order=3)

# Create a discrete colormap
cmap = ListedColormap(colors[:5])

# Define boundaries between color bins
bounds = np.arange(-1.5, 4.5, 1)  # [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
norm = BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(8, 6))

# pcm = ax.pcolormesh(X, Y, PFEI_min_nacelle_array, cmap='viridis')
ax.imshow(
    smooth_field, extent=[X.min() / 1e6, X.max() / 1e6, Y.min() / 1e3, Y.max() / 1e3],
    cmap=cmap, norm=norm, origin='lower', aspect='auto',
)

# ax.contourf(X, Y, min_idx, levels=[0, 1, 2, 3], cmap=cmap, norm=norm)

ax.set_xlabel('FCS power density (kW/l)', labelpad=10)
ax.set_ylabel('FCS specific power (kW/kg)', labelpad=10)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

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

#%%

fig, ax = plt.subplots()

ax.scatter(P_prop_max_nacelle_list[0] / Vol_nacelle_list[0] / 0.55 / 1e6, PFEI_nacelle_list[0], alpha=0.1)
ax.scatter(P_prop_max_wing_list[0] / Vol_wing_list[0] / 0.36 / 1e6, PFEI_wing_list[0], alpha=0.1)
ax.scatter(P_prop_max_fuselage_rear_list[0] / Vol_fuselage_rear_list[0] / 0.69 / 1e6, PFEI_fuselage_rear_list[0], alpha=0.1)
ax.scatter(P_prop_max_fuselage_underfloor_list[0] / Vol_fuselage_underfloor_list[0] / 0.88 / 1e6, PFEI_fuselage_underfloor_list[0], alpha=0.1, zorder=-100)

ax.scatter(P_prop_max_nacelle_list[1] / Vol_nacelle_list[1] / 0.55 / 1e6, PFEI_nacelle_list[1], alpha=0.1)
ax.scatter(P_prop_max_wing_list[1] / Vol_wing_list[1] / 0.36 / 1e6, PFEI_wing_list[1], alpha=0.1)
ax.scatter(P_prop_max_fuselage_rear_list[1] / Vol_fuselage_rear_list[1] / 0.69 / 1e6, PFEI_fuselage_rear_list[1], alpha=0.1)
ax.scatter(P_prop_max_fuselage_underfloor_list[1] / Vol_fuselage_underfloor_list[1] / 0.88 / 1e6, PFEI_fuselage_underfloor_list[1], alpha=0.1, zorder=-100)

ax.set_xlim(0, 10)

plt.show()


