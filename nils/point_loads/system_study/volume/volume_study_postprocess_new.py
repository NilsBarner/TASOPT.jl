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

# # xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'volume_results_narrowbody_291025_v1.xlsx'))
# # xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'volume_results_narrowbody_fuselage.xlsx'))
# # xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'volume_results_narrowbody_wing.xlsx'))
# # xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'volume_results_narrowbody_nacelle.xlsx'))
# xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'volume_results_narrowbody_all.xlsx'))

# sheet_dict = {}
# for sheet_name in xlsx.sheet_names:
#     sheet_dict[sheet_name] = pd.read_excel(xlsx, sheet_name=sheet_name)
    
# df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_wing_newest.csv'))
df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_nacelle_newest.csv'))
# df = pd.read_csv(os.path.join(os.getcwd(), 'volume_results_narrowbody_fuselage_newest.csv'))

# sys.exit()

#%%

# df = sheet_dict['results']

fcs_loc_unique = df["fcs_loc"].unique()
radius_unique = np.sort(df["radius"].unique())
AR_unique = np.sort(df["AR"].unique())
tdivc_scale_unique = np.sort(df["tdivc_scale"].unique())
N_eng_unique = np.sort(df["N_eng"].unique())
HTR_f_unique = np.sort(df["HTR_f"].unique())
Vspec_unique = np.sort(df["Vspec"].unique())
fcs_fuselage_location_unique = np.sort(df["fcs_fuselage_location"].unique())

# Create 2D coordinate grids
fcs_loc_grid, radius_grid, AR_grid, tdivc_scale_grid, N_eng_grid, HTR_f_grid, Vspec_grid, fcs_fuselage_location_grid = np.meshgrid(
    fcs_loc_unique, radius_unique, AR_unique, tdivc_scale_unique, N_eng_unique, HTR_f_unique, Vspec_unique, fcs_fuselage_location_unique, indexing="ij"
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

# Initialize empty arrays
CDS_grid = np.full(shape, np.nan)
CDS_ref_grid = np.full(shape, np.nan)
WMTO_grid = np.full(shape, np.nan)
WMTO_ref_grid = np.full(shape, np.nan)
Vol_wing_grid = np.full(shape, np.nan)
PFEI_grid = np.full(shape, np.nan)
PFEI_ref_grid = np.full(shape, np.nan)
seats_abreast_grid = np.full(shape, np.nan)
seats_abreast_ref_grid = np.full(shape, np.nan)
Vol_nacelle_grid = np.full(shape, np.nan)
y_centroid_wing_grid = np.full(shape, np.nan)
L_fuse_grid = np.full(shape, np.nan)
y_centroid_nacelles_grid = np.full(shape, np.nan)
span_grid = np.full(shape, np.nan)

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
    
    CDS_grid[h, i, j, k, l, m, n, o] = row["CDS"]
    # CDS_ref_grid[h, i, j, k, l, m, n, o] = row["CDS_ref"]
    WMTO_grid[h, i, j, k, l, m, n, o] = row["WMTO"]
    # WMTO_ref_grid[h, i, j, k, l, m, n, o] = row["WMTO_ref"]
    Vol_wing_grid[h, i, j, k, l, m, n, o] = row["Vol_wing"]
    PFEI_grid[h, i, j, k, l, m, n, o] = row["PFEI"]
    # PFEI_ref_grid[h, i, j, k, l, m, n, o] = row["PFEI_ref"]
    seats_abreast_grid[h, i, j, k, l, m, n, o] = row["seats_abreast"]
    # seats_abreast_ref_grid[h, i, j, k, l, m, n, o] = row["seats_abreast_ref"]
    Vol_nacelle_grid[h, i, j, k, l, m, n, o] = row["Vol_nacelle"]
    y_centroid_wing_grid[h, i, j, k, l, m, n, o] = row["y_centroid_wing"]
    # L_fuse_grid[h, i, j, k, l, m, n, o] = row["L_fuse"]
    y_centroid_nacelles_grid[h, i, j, k, l, m, n, o] = row["y_centroid_nacelles"]
    span_grid[h, i, j, k, l, m, n, o] = row["span"]
    
# sys.exit()

y_centroid_wing_norm_grid = y_centroid_wing_grid / (span_grid / 2)
y_centroid_nacelles_norm_grid = y_centroid_nacelles_grid / (span_grid / 2)

#%% Wing
"""
mask = (0, slice(None), slice(None), slice(None), slice(None), slice(None), 0, slice(None))

fig, ax = plt.subplots(figsize=tuple(np.array([6.4, 4.8]) * 1.5))

# s = ax.scatter(
#     WMTO_grid[mask],# / WMTO_ref_grid[mask],
#     # CDS_grid[mask] / CDS_ref_grid[mask],
#     PFEI_grid[mask],# / PFEI_ref_grid[mask],
#     # c=radius_grid[mask],
#     # c=N_eng_grid[mask],
#     # c=HTR_f_grid[mask],
#     # c=Vspec_grid[mask],
#     # c=Vol_wing_grid[mask],
#     # c=AR_grid[mask],
#     # c=tdivc_scale_grid[mask],
#     marker='.',
#     color='black',
#     zorder=100,
# )

# ax.contour(
#     np.squeeze(WMTO_grid[mask]),
#     np.squeeze(PFEI_grid[mask]),
#     np.squeeze(AR_grid[mask]),
#     levels=AR_unique,
#     colors=['blue'],
# )

# ax.contour(
#     np.squeeze(WMTO_grid[mask]),
#     np.squeeze(PFEI_grid[mask]),
#     np.squeeze(tdivc_scale_grid[mask]),
#     levels=tdivc_scale_unique,
#     colors=['red'],
# )
'''
import alphashape
from alpha_shape_creator import ensure_xy_array, add_shapely_polygon

finite_mask = (np.isfinite(WMTO_grid[mask]) & np.isfinite(PFEI_grid[mask]))

points_2d = np.array([
    WMTO_grid[mask][finite_mask].flatten(),
    PFEI_grid[mask][finite_mask].flatten(),
]).T

pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
alpha_shape = alphashape.alphashape(pts, 10)
add_shapely_polygon(
    ax, alpha_shape, facecolor=colors[0], edgecolor='None', alpha=0.2,
)

# =============================================================================
cvhull_points = points_2d
if len(cvhull_points) > 0:
    cvhull = ConvexHull(cvhull_points)
    cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
    polygon = Polygon(cvhull_vertices)
    x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    ax.fill(
        x_outline_polygon, y_outline_polygon,
        hatch='..', facecolor=colors[0], edgecolor='None', alpha=0.2,
    )
# =============================================================================
'''
for i, AR in enumerate(AR_unique):
    
    x_slice = WMTO_grid[0, 0, i, slice(None), 0, 0, 0, 0]
    y_slice = PFEI_grid[0, 0, i, slice(None), 0, 0, 0, 0]
    
    good = np.isfinite(x_slice) & np.isfinite(y_slice)
    
    ax.plot(
        x_slice[good], y_slice[good],
        color=colors[0],
    )
        
for j, tdivc_scale in enumerate(tdivc_scale_unique):
    
    good = np.isfinite(x_slice) & np.isfinite(y_slice)
    
    x_slice = WMTO_grid[0, 0, slice(None), j, 0, 0, 0, 0]
    y_slice = PFEI_grid[0, 0, slice(None), j, 0, 0, 0, 0]
    
    good = np.isfinite(x_slice) & np.isfinite(y_slice)
    
    ax.plot(
        x_slice[good], y_slice[good],
        color=colors[1],
    )

ax.set_xlabel(r'Drag', labelpad=10)
ax.set_ylabel(r'Weight', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_elements = [
    Line2D([0], [0], color=colors[0], alpha=1.0, label=r'Contours of aspect ratio'),
    Line2D([0], [0], color=colors[1], alpha=1.0, label=r'Contours of thickness-chord ratio'),
]
leg_main = ax.legend(handles=legend_elements, loc='upper left', frameon=False)
ax.add_artist(leg_main)

add_margin(ax, m=0.05)

# fig.colorbar(s, cmap='viridis')

plt.show()

#%%

P_req = 20e6
Vspec_curr = 0.85e6
Vol_wing_req = P_req / Vspec_curr

Vol_wing_norm_grid = Vol_wing_grid / Vol_wing_req

fig, ax = plt.subplots(figsize=tuple(np.array([6.4, 4.8]) * 1.5))

# s = ax.scatter(
#     Vol_wing_norm_grid[mask],# / WMTO_ref_grid[mask],
#     # CDS_grid[mask] / CDS_ref_grid[mask],
#     y_centroid_wing_norm_grid[mask],# / PFEI_ref_grid[mask],
#     # c=radius_grid[mask],
#     # c=N_eng_grid[mask],
#     # c=HTR_f_grid[mask],
#     # c=Vspec_grid[mask],
#     # c=Vol_wing_norm_grid[mask],
#     # c=AR_grid[mask],
#     # c=tdivc_scale_grid[mask],
#     marker='.',
#     color='black',
#     zorder=100,
# )

# ax.contour(
#     np.squeeze(Vol_wing_norm_grid[mask]),
#     np.squeeze(y_centroid_wing_norm_grid[mask]),
#     np.squeeze(AR_grid[mask]),
#     levels=AR_unique,
#     colors=['blue'],
# )

# ax.contour(
#     np.squeeze(Vol_wing_norm_grid[mask]),
#     np.squeeze(y_centroid_wing_norm_grid[mask]),
#     np.squeeze(tdivc_scale_grid[mask]),
#     levels=tdivc_scale_unique,
#     colors=['red'],
# )

import alphashape
from alpha_shape_creator import ensure_xy_array, add_shapely_polygon

finite_mask = (np.isfinite(Vol_wing_norm_grid[mask]) & np.isfinite(y_centroid_wing_norm_grid[mask]))

points_2d = np.array([
    Vol_wing_norm_grid[mask][finite_mask].flatten(),
    y_centroid_wing_norm_grid[mask][finite_mask].flatten(),
]).T

# pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
# alpha_shape = alphashape.alphashape(pts, 0.1)
# add_shapely_polygon(
#     ax, alpha_shape, facecolor=colors[0], edgecolor='None', alpha=0.2,
# )

# =============================================================================
cvhull_points = points_2d
if len(cvhull_points) > 0:
    cvhull = ConvexHull(cvhull_points)
    cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
    polygon = Polygon(cvhull_vertices)
    x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    ax.fill(
        x_outline_polygon, y_outline_polygon,
        hatch='..', facecolor=colors[2], edgecolor='None', alpha=0.3,
    )
# =============================================================================

for i, AR in enumerate(AR_unique):
    
    x_slice = Vol_wing_norm_grid[0, 0, i, slice(None), 0, 0, 0, 0]
    y_slice = y_centroid_wing_norm_grid[0, 0, i, slice(None), 0, 0, 0, 0]
    
    good = np.isfinite(x_slice) & np.isfinite(y_slice)
    
    ax.plot(
        x_slice[good], y_slice[good],
        color=colors[0],
    )
        
for j, tdivc_scale in enumerate(tdivc_scale_unique):
    
    good = np.isfinite(x_slice) & np.isfinite(y_slice)
    
    x_slice = Vol_wing_norm_grid[0, 0, slice(None), j, 0, 0, 0, 0]
    y_slice = y_centroid_wing_norm_grid[0, 0, slice(None), j, 0, 0, 0, 0]
    
    good = np.isfinite(x_slice) & np.isfinite(y_slice)
    
    ax.plot(
        x_slice[good], y_slice[good],
        color=colors[1],
    )

# fig.colorbar(s, cmap='viridis')

ax.set_xlabel(r'Available / required FCS volume (-)', labelpad=10)
ax.set_ylabel(r'Span-normalised location of centroid (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_elements = [
    Line2D([0], [0], color=colors[0], alpha=1.0, label=r'Aspect ratio'),
    Line2D([0], [0], color=colors[1], alpha=1.0, label=r'Thickness-chord ratio'),
    Patch(facecolor=colors[2], edgecolor='None', alpha=0.3, label=r'Insufficient FCS volume'),
]
leg_main = ax.legend(handles=legend_elements, loc='center right', frameon=False)
ax.add_artist(leg_main)

add_margin(ax, m=0.05)

# plt.savefig('wing_centroid_vs_volume.svg', format='svg')
plt.show()

sys.exit()
"""
#%% Nacelle
# '''

P_req = 20e6
Vspec_curr = 0.85e6
Vol_nacelle_req = P_req / Vspec_curr

Vol_nacelle_norm_grid = Vol_nacelle_grid / Vol_nacelle_req

mask = (0, 0, 0, 0, slice(None), slice(None), 0, 0)

fig, ax = plt.subplots(figsize=tuple(np.array([6.4, 4.8]) * 1.5))

s = ax.scatter(
    Vol_nacelle_norm_grid[0, 0, 0, 0, 0, slice(None), 0, 0],
    # CDS_grid[mask] / CDS_ref_grid[mask],
    y_centroid_nacelles_norm_grid[0, 0, 0, 0, 0, slice(None), 0, 0],
    # c=N_eng_grid[mask],
    # c=HTR_f_grid[mask],
    marker='.',# color='blue'
    color=colors[1],
    zorder=100,
)

# s = ax.contour(
#     Vol_nacelle_norm_grid[mask],
#     y_centroid_nacelles_norm_grid[mask],
#     N_eng_grid[mask],
#     levels=N_eng_unique,
#     colors=['blue'],
# )

# s = ax.contour(
#     Vol_nacelle_norm_grid[mask],
#     y_centroid_nacelles_norm_grid[mask],
#     HTR_f_grid[mask],
#     levels=HTR_f_unique,
#     colors=['red'],
# )

cvhull_points = np.array([
    Vol_nacelle_norm_grid[0, 0, 0, 0, slice(None), slice(None), 0, 0].flatten(),
    y_centroid_nacelles_norm_grid[0, 0, 0, 0, slice(None), slice(None), 0, 0].flatten(),
]).T
if len(cvhull_points) > 0:
    cvhull = ConvexHull(cvhull_points)
    cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
    polygon = Polygon(cvhull_vertices)
    x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    ax.fill(
        x_outline_polygon, y_outline_polygon,
        hatch='..', facecolor=colors[2], edgecolor='None', alpha=0.3,
    )

for i in range(len(N_eng_unique)):
    ax.plot(
        Vol_nacelle_norm_grid[0, 0, 0, 0, i, slice(None), 0, 0],
        y_centroid_nacelles_norm_grid[0, 0, 0, 0, i, slice(None), 0, 0],
        color=colors[0],
    )
    print(y_centroid_nacelles_norm_grid[0, 0, 0, 0, i, slice(None), 0, 0][-1])
    
for j in range(len(HTR_f_unique)):
    ax.plot(
        Vol_nacelle_norm_grid[0, 0, 0, 0, slice(None), j, 0, 0][1:],
        y_centroid_nacelles_norm_grid[0, 0, 0, 0, slice(None), j, 0, 0][1:],
        color=colors[1],
    )

ax.set_xlabel(r'Available / required FCS volume (-)', labelpad=10)
ax.set_ylabel(r'Span-normalised location of centroid (-)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_elements = [
    Line2D([0], [0], color=colors[0], alpha=1.0, label=r'Number of engines'),
    Line2D([0], [0], color=colors[1], alpha=1.0, label=r'Fan tub-tip-ratio'),
    Patch(facecolor=colors[2], edgecolor='None', alpha=0.3, label=r'Insufficient FCS volume'),
]
leg_main = ax.legend(handles=legend_elements, loc='center right', frameon=False)
ax.add_artist(leg_main)

add_margin(ax, m=0.05)

# fig.colorbar(s, cmap='viridis')

# plt.savefig('nacelle_centroid_vs_volume.svg', format='svg')
plt.show()

#%%

mask = (0, 0, 0, 0, slice(None), slice(None), 0, 0)

fig, ax = plt.subplots()

# s = ax.scatter(
#     WMTO_grid[mask],
#     CDS_grid[mask],
#     # c=N_eng_grid[mask],
#     # c=HTR_f_grid[mask],
#     marker='.',# color='blue'
# )

s = ax.contour(
    WMTO_grid[mask],
    CDS_grid[mask],
    N_eng_grid[mask],
    levels=N_eng_unique,
    colors=['blue'],
)

s = ax.contour(
    WMTO_grid[mask],
    CDS_grid[mask],
    HTR_f_grid[mask],
    levels=HTR_f_unique,
    colors=['red'],
)

fig.colorbar(s, cmap='viridis')

plt.show()

sys.exit()
# '''
#%% Fuselage

from matplotlib.colors import Normalize
from matplotlib_custom_settings import *
import matplotlib.colors as mcolors

def opaque_color_from_hex(hex_color, alpha=0.5, background='white'):
    """Return an opaque RGB color that looks like `hex_color` with transparency `alpha`
    over the given `background`."""
    rgba_fg = np.array(mcolors.to_rgba(hex_color))
    rgba_bg = np.array(mcolors.to_rgba(background))
    blended_rgb = rgba_fg[:3] * alpha + rgba_bg[:3] * (1 - alpha)
    return mcolors.to_hex(blended_rgb, keep_alpha=False)


# your masks
mask_0 = (0, slice(None), 0, 0, 0, 0, slice(None), 0)
mask_1 = (0, slice(None), 0, 0, 0, 0, slice(None), 1)

# Example: assuming seats_abreast_grid is your color variable
# vmin = np.nanmin(seats_abreast_grid)
# vmax = np.nanmax(seats_abreast_grid)
# vmin = np.nanmin(Vspec_grid)
# vmax = np.nanmax(Vspec_grid)
vmin = np.nanmin(radius_grid)
vmax = np.nanmax(radius_grid)
norm = Normalize(vmin=vmin, vmax=vmax)

fig, ax = plt.subplots(figsize=tuple(np.array([6.4, 4.8]) * 1.5))

# # optional: plot outlines or reference markers
# ax.scatter(radius_grid[mask_1], L_fuse_grid[mask_1], color='red')
# ax.scatter(radius_grid[mask_0], L_fuse_grid[mask_0], color='blue')

# fig2, ax2 = plt.subplots()

for i, seats_abreast in enumerate(np.unique(seats_abreast_grid)):
    
    mask_seats_abreast = (seats_abreast_grid == seats_abreast) & (fcs_fuselage_location_grid == 1.0)
    
    radius_grid_masked = radius_grid[mask_seats_abreast]
    L_fuse_grid_masked = L_fuse_grid[mask_seats_abreast]
    finite_mask = np.isfinite(radius_grid_masked) & np.isfinite(L_fuse_grid_masked)
    radius_grid_filtered_masked = radius_grid_masked[finite_mask]
    L_fuse_grid_filtered_masked = L_fuse_grid_masked[finite_mask]
    
    cvhull_points = np.array([
        radius_grid_filtered_masked.flatten(),
        L_fuse_grid_filtered_masked.flatten(),
    ]).T
    if len(cvhull_points) > 0:
        cvhull = ConvexHull(cvhull_points)
        cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
        polygon = Polygon(cvhull_vertices)
        x_outline_polygon, y_outline_polygon = polygon.exterior.xy
        ax.fill(
            x_outline_polygon, y_outline_polygon,
            hatch='..', facecolor=opaque_color_from_hex(colors[i], 0.5), edgecolor=opaque_color_from_hex(colors[i], 1), zorder=-1,
        )
        
    ###
        
    x_masked = np.ma.array(radius_grid, mask=~mask_seats_abreast)
    y_masked = np.ma.array(L_fuse_grid, mask=~mask_seats_abreast)
    z_masked = np.ma.array(Vspec_grid, mask=~mask_seats_abreast)
        
    # ax.contour(
    #     x_masked[mask_1],
    #     y_masked[mask_1],
    #     z_masked[mask_1],
    # )
    
    # for j, Vspec in enumerate(Vspec_unique):
    #     ax.plot(
    #         x_masked[0, slice(None), 0, 0, 0, 0, j, 1],
    #         y_masked[0, slice(None), 0, 0, 0, 0, j, 1],
    #         color='black',
    #     )
    
    for j, Vspec in enumerate(Vspec_unique):
        
        if not (j == 0 or j == len(Vspec_unique) - 1):
            continue
        
        # take the 1-D slice (masked array) you already built
        x_slice = x_masked[0, slice(None), 0, 0, 0, 0, j, 1]
        y_slice = y_masked[0, slice(None), 0, 0, 0, 0, j, 1]
    
        # get boolean mask arrays (True where masked)
        x_mask_bool = np.ma.getmaskarray(x_slice)
        y_mask_bool = np.ma.getmaskarray(y_slice)
    
        # construct boolean of good points: not masked AND finite
        x_data = x_slice.data  # view of underlying data
        y_data = y_slice.data
        good = (~x_mask_bool) & (~y_mask_bool) & np.isfinite(x_data) & np.isfinite(y_data)
    
        # if no good points, skip
        if not np.any(good):
            continue
    
        # select only good entries (preserves order so neighbors connect)
        x_plot = x_data[good]
        y_plot = y_data[good]
    
        # =============================================================================
        ax.plot(x_plot, y_plot, color='black')
        # =============================================================================
        
        # L_fuse_0 = L_fuse_grid[0, slice(None), 0, 0, 0, 0, j, 0]
        # L_fuse_1 = L_fuse_grid[0, slice(None), 0, 0, 0, 0, j, 1]
        
        # gt_mask = L_fuse_0 > L_fuse_1
        # ax.scatter(x_slice[gt_mask], y_slice[gt_mask], color='k', marker='.')
        
        # print('x_slice[gt_mask], y_slice[gt_mask] =', x_slice[gt_mask], y_slice[gt_mask])
        
        # try:
        #     cvhull_points = np.array([
        #         x_slice[gt_mask],
        #         y_slice[gt_mask],
        #     ]).T
        #     cvhull = ConvexHull(cvhull_points)
        #     cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
        #     polygon = Polygon(cvhull_vertices)
        #     x_outline_polygon, y_outline_polygon = polygon.exterior.xy
        #     ax.fill(
        #         x_outline_polygon, y_outline_polygon, color=colors[i], alpha=0.5, zorder=-1,
        #     )
        # except Exception as e:
        #     # print(str(e))
        #     pass
        
    # =============================================================================
    # L_fuse_0 = L_fuse_grid[0, slice(None), 0, 0, 0, 0, slice(None), 0]
    # L_fuse_1 = L_fuse_grid[0, slice(None), 0, 0, 0, 0, slice(None), 1]
    
    # gt_mask = L_fuse_0 > L_fuse_1
    # ax.scatter(x_masked[mask_1][gt_mask], y_masked[mask_1][gt_mask], color=colors[i], marker='.')
    
    # print(x_masked[mask_1][gt_mask], y_masked[mask_1][gt_mask])
    # print()
    
    
    _mask_seats_abreast = (seats_abreast_grid == seats_abreast) & (fcs_fuselage_location_grid == 0.0)
    _L_fuse_grid_masked = L_fuse_grid[_mask_seats_abreast]
    # _finite_mask = np.isfinite(radius_grid_masked) & np.isfinite(L_fuse_grid_masked)
    # print('seats_abreast =', seats_abreast)
    
    
    try:
        _gt_mask = _L_fuse_grid_masked > L_fuse_grid_masked
        
        # ax.scatter(
        #     radius_grid_filtered_masked[_gt_mask].data,
        #     L_fuse_grid_filtered_masked[_gt_mask].data,
        #     color=colors[i], marker='.',
        # )
        
        # # =============================================================================
        # Vspec_grid_masked = Vspec_grid[mask_seats_abreast]
        # finite_mask = np.isfinite(radius_grid_masked) & np.isfinite(Vspec_grid_masked)
        # Vspec_grid_filtered_masked = Vspec_grid_masked[finite_mask]
        
        # ax2.scatter(
        #     radius_grid_filtered_masked[_gt_mask].data,
        #     Vspec_grid_filtered_masked[_gt_mask].data,
        #     color=colors[i], marker='.',
        # )
        # # =============================================================================
        
        cvhull_points = np.array([
            # x_masked[mask_1][gt_mask].data,
            # y_masked[mask_1][gt_mask].data,
            radius_grid_filtered_masked[_gt_mask].data,
            L_fuse_grid_filtered_masked[_gt_mask].data,
        ]).T
        cvhull = ConvexHull(cvhull_points)
        cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
        polygon = Polygon(cvhull_vertices)
        x_outline_polygon, y_outline_polygon = polygon.exterior.xy
        ax.fill(
            x_outline_polygon,
            y_outline_polygon,
            hatch='\\\\', facecolor=opaque_color_from_hex(colors[i], 0.5), edgecolor=opaque_color_from_hex(colors[i], 1), zorder=-1,
        )
    except Exception as e:
        # print(str(e))
        pass
    
        # try:
        #     _gt_mask = _L_fuse_grid_masked > L_fuse_grid_masked
        #     ax.scatter(
        #         radius_grid_filtered_masked[_gt_mask].data,
        #         L_fuse_grid_filtered_masked[_gt_mask].data,
        #         color=colors[i], marker='.',
        #     )
        # except:
        #     pass
    # =============================================================================

        
    # if seats_abreast == 7:
    #     print(x_masked[0, slice(None), 0, 0, 0, 0, slice(None), 1])
    #     sys.exit('Stop.')
    
    ###
    
    mask_seats_abreast = ((seats_abreast_grid == seats_abreast) & (fcs_fuselage_location_grid == 0.0))
    try:
        L_fuses_binary = np.unique(L_fuse_grid[mask_seats_abreast])
        print('L_fuses_binary =', L_fuses_binary)
        mask_top = (L_fuse_grid[mask_seats_abreast] == np.nanmax(L_fuses_binary))
        mask_bottom = (L_fuse_grid[mask_seats_abreast] == np.nanmin(L_fuses_binary))
        ax.plot(
            radius_grid[mask_seats_abreast][mask_top],
            L_fuse_grid[mask_seats_abreast][mask_top],
            color=colors[i],
            marker='s',
        )
        
        ax.plot(
            radius_grid[mask_seats_abreast][mask_bottom],
            L_fuse_grid[mask_seats_abreast][mask_bottom],
            color=colors[i],
            marker='s',
        )
    except:
        pass
    
    # radius_grid_masked = radius_grid[mask_seats_abreast]
    # L_fuse_grid_masked = L_fuse_grid[mask_seats_abreast]
    # finite_mask = np.isfinite(radius_grid_masked) & np.isfinite(L_fuse_grid_masked)
    # radius_grid_filtered_masked = radius_grid_masked[finite_mask]
    # L_fuse_grid_filtered_masked = L_fuse_grid_masked[finite_mask]
    
    # cvhull_points = np.array([
    #     radius_grid_filtered_masked.flatten(),
    #     L_fuse_grid_filtered_masked.flatten(),
    # ]).T
    # try:
    #     cvhull = ConvexHull(cvhull_points)
    #     cvhull_vertices = cvhull_points[cvhull.vertices]                  # ordered hull vertices
    #     polygon = Polygon(cvhull_vertices)
    #     x_outline_polygon, y_outline_polygon = polygon.exterior.xy
    #     ax.fill(
    #         x_outline_polygon, y_outline_polygon, color=colors[i], alpha=0.5, zorder=-1,
    #         linestyle='dashed',
    #     )
    # except Exception as e:
    #     # print(str(e))
    #     pass
    
    
    

# # these will share the same color scale
# s1 = ax.scatter(
#     radius_grid[mask_1],
#     L_fuse_grid[mask_1],
#     # c=seats_abreast_grid[mask_1],
#     # c=Vspec_grid[mask_1],
#     c=radius_grid[mask_1],
#     marker='.',
#     norm=norm,
#     cmap='viridis'
# )
# s2 = ax.scatter(
#     radius_grid[mask_0],
#     L_fuse_grid[mask_0],
#     # c=seats_abreast_grid[mask_0],
#     # c=Vspec_grid[mask_0],
#     c=radius_grid[mask_0],
#     marker='o',
#     norm=norm,
#     cmap='viridis'
# )

# # only one colorbar is needed — use the last scatter as representative
# fig.colorbar(s2, ax=ax)

ax.set_xlabel('Fuselage radius (m)', labelpad=10)
ax.set_ylabel('Fuselage length (m)', labelpad=10)
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

legend_elements = [
    Line2D([0], [0], color='black', marker='s', alpha=1.0, label=r'Underfloor FCS'),
    Patch(facecolor='None', edgecolor='black', alpha=1.0, label='Aft FCS'),
    Patch(facecolor='None', edgecolor='black', hatch='..', alpha=1.0, label=r'$L_\mathrm{fuse,uf} < L_\mathrm{fuse,aft}$ for same $\nu_\mathrm{FCS}$'),
    Patch(facecolor='None', edgecolor='black', hatch='\\\\', alpha=1.0, label=r'$L_\mathrm{fuse,uf} > L_\mathrm{fuse,aft}$ for same $\nu_\mathrm{FCS}$'),
]
leg_main = ax.legend(handles=legend_elements, loc='upper right', frameon=False)
ax.add_artist(leg_main)

# plt.savefig('fuse_length_vs_radius.svg', format='svg')
plt.show()



