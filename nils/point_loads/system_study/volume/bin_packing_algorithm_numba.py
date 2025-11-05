#!/usr/bin/env python3
"""
wing_box_packer_numba.py

Robust, physically-consistent cuboid packer for a wing-box volume with:
 - strict full-containment checks (cuboid's 8 corners inside mesh)
 - optional 90° rotations about Z (swap X/Y)
 - Numba-accelerated greedy placement (no Python fallback)

Requirements:
    pip install numpy scipy trimesh matplotlib numba

Author: ChatGPT (corrected & numba-accelerated)
"""
import numpy as np
from scipy.signal import convolve
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numba import njit

# -------------------------
# Geometry: parametric wing-box (unchanged)
# -------------------------
def make_tapered_wing_box(span=10.0, root_chord=3.0, tip_chord=1.2,
                          thickness=0.3, sweep_le=0.0, dihedral_deg=0.0):
    half = span / 2.0
    ys = np.array([-half, half])
    chords = np.array([root_chord, tip_chord])
    le_x = np.array([0.0, sweep_le])
    dihedral_rad = np.deg2rad(dihedral_deg)
    zs = np.array([0.0, np.tan(dihedral_rad) * half])

    verts = []
    for j in range(2):
        y = ys[j]
        c = chords[j]
        le = le_x[j]
        z_off = zs[j]
        verts.append([le,       y, z_off + thickness/2.0])   # top leading
        verts.append([le + c,   y, z_off + thickness/2.0])   # top trailing
        verts.append([le + c,   y, z_off - thickness/2.0])   # bottom trailing
        verts.append([le,       y, z_off - thickness/2.0])   # bottom leading

    def v(i,j): return 4*i + j
    faces = []
    faces += [[v(0,0), v(1,0), v(1,1)], [v(0,0), v(1,1), v(0,1)]]  # top
    faces += [[v(0,3), v(1,3), v(1,2)], [v(0,3), v(1,2), v(0,2)]]  # bottom
    faces += [[v(0,0), v(0,3), v(1,3)], [v(0,0), v(1,3), v(1,0)]]  # leading edge
    faces += [[v(0,1), v(1,1), v(1,2)], [v(0,1), v(1,2), v(0,2)]]  # trailing

    verts = np.array(verts)
    mesh = trimesh.Trimesh(vertices=verts, faces=np.array(faces), process=True)

    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        mesh = trimesh.Trimesh(vertices=verts).convex_hull
    return mesh

# =============================================================================
# def create_mesh_from_coordinates(coords):
#     xs = coords['xw']
#     ys = coords['yw']
#     zs = coords['zw']
    
#     verts = np.vstack((xs, ys, zs)).T
#     print('verts =', verts)
    
#     faces = []
#     faces += [verts[0], verts[1], verts[6], verts[7]]
#     faces += [verts[1], verts[2], verts[7], verts[8]]
    
#     faces += [verts[2], verts[3], verts[8], verts[9]]
    
#     faces += [verts[6], verts[7], verts[10], verts[11]]
#     faces += [verts[7], verts[8], verts[9], verts[10]]
    
#     faces += [verts[4], verts[5], verts[10], verts[11]]
#     faces += [verts[3], verts[4], verts[9], verts[10]]
    
#     faces += [verts[0], verts[1], verts[4], verts[5]]
#     faces += [verts[1], verts[2], verts[3], verts[4]]
    
#     mesh = trimesh.Trimesh(
#         vertices=verts, faces=np.array(faces), process=True
#     )

#     if not mesh.is_watertight:
#         try:
#             mesh.remove_unreferenced_vertices()
#             mesh.fill_holes()
#         except Exception:
#             pass
#     if not mesh.is_watertight:
#         mesh = trimesh.Trimesh(vertices=verts).convex_hull
#     return mesh

import numpy as np
import trimesh

def create_mesh_from_coordinates(coords):
    xs = coords['xw'].to_numpy()
    ys = coords['yw'].to_numpy()
    zs = coords['zw'].to_numpy()

    verts = np.vstack((xs, ys, zs)).T

    # define quads by vertex indices (not coordinates)
    quads = [
        [0, 1, 6, 7],
        [1, 2, 7, 8],
        [2, 3, 8, 9],
        [6, 7, 10, 11],
        [7, 8, 9, 10],
        [4, 5, 10, 11],
        [3, 4, 9, 10],
        [0, 1, 4, 5],
        [1, 2, 3, 4],
    ]

    # convert each quad to two triangles
    faces = []
    for a, b, c, d in quads:
        faces.append([a, b, c])
        faces.append([a, c, d])

    faces = np.array(faces, dtype=np.int64)

    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=True)

    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        mesh = trimesh.Trimesh(vertices=verts).convex_hull

    return mesh

# =============================================================================

# -------------------------
# Numba-accelerated greedy placement
# -------------------------
@njit
def greedy_place_numba(occ, nx, ny, nz,
                       xi0, xi1, yi0, yi1, zi0, zi1,
                       order_idx):
    """
    occ : boolean 3D array (nx,ny,nz) modified in-place
    xi0/xi1 etc : int arrays (length M) of voxel-index ranges [xi0, xi1) ...
    order_idx : int array (length M) specifying the order in which to try candidates
    Returns: placed_flags boolean array of length M indicating which candidates were placed
    """
    M = order_idx.shape[0]
    placed = np.zeros(M, dtype=np.bool_)
    for p in range(M):
        cand = order_idx[p]
        x0 = xi0[cand]; x1 = xi1[cand]
        y0 = yi0[cand]; y1 = yi1[cand]
        z0 = zi0[cand]; z1 = zi1[cand]

        # quick reject for empty ranges
        if x1 <= x0 or y1 <= y0 or z1 <= z0:
            continue

        overlap = False
        for i in range(x0, x1):
            for j in range(y0, y1):
                for k in range(z0, z1):
                    if occ[i, j, k]:
                        overlap = True
                        break
                if overlap:
                    break
            if overlap:
                break
        if not overlap:
            # place it (mark voxels occupied)
            placed[cand] = True
            for i in range(x0, x1):
                for j in range(y0, y1):
                    for k in range(z0, z1):
                        occ[i, j, k] = True
    return placed

# -------------------------
# Core packer (robust + rotation)
# -------------------------
def pack_cuboids_in_mesh(mesh, cuboid_dims, pitch=None, allow_z90=False, max_voxels=200**3, pitch_tol=1e-8):
    """
    Pack axis-aligned cuboids into mesh interior.
    - cuboid_dims: (Lx, Ly, Lz) in meters (Lx along X/flight direction).
    - pitch: voxel pitch. If None, a safe automatic pitch is chosen.
      If user supplies pitch, it must satisfy either:
        * pitch <= min_dim / 2      (high resolution)  OR
        * each (dim / pitch) is within pitch_tol of an integer (exact multiple)
      Otherwise ValueError is raised to avoid incorrect placements.
    - allow_z90: if True, also allow cuboid rotated by 90 deg about Z (swap X/Y).
    Returns dict similar to earlier versions; each placement includes orientation (0 or 1).
    """
    Lx, Ly, Lz = cuboid_dims
    if min(Lx, Ly, Lz) <= 0:
        raise ValueError("cuboid_dims must be positive")

    min_dim = min(Lx, Ly, Lz)

    # Choose or validate pitch
    if pitch is None:
        # heuristic: at least ~3 voxels on smallest side, but ensure voxels not huge
        pitch = min_dim / 3.0
    else:
        pitch = float(pitch)

    # If user provided pitch, validate it to avoid pathological cases.
    # Either pitch is small enough (<= min_dim/2) OR dims are integer multiples of pitch (within tol).
    def is_integer_multiple(val, p):
        return abs(round(val / p) - (val / p)) <= pitch_tol

    if not (pitch <= (min_dim / 2.0) or (is_integer_multiple(Lx, pitch) and is_integer_multiple(Ly, pitch) and is_integer_multiple(Lz, pitch))):
        raise ValueError(
            "Provided pitch is incompatible with cuboid dimensions. "
            "Either choose a smaller pitch (<= min_dim/2) or a pitch that divides each cuboid dimension exactly."
        )

    # Build voxel centers from mesh bounding box
    bb_min = mesh.bounds[0].astype(float)
    bb_max = mesh.bounds[1].astype(float)

    xs = np.arange(bb_min[0] + 0.5 * pitch, bb_max[0], pitch)
    ys = np.arange(bb_min[1] + 0.5 * pitch, bb_max[1], pitch)
    zs = np.arange(bb_min[2] + 0.5 * pitch, bb_max[2], pitch)

    nx, ny, nz = len(xs), len(ys), len(zs)
    total_vox = int(nx) * int(ny) * int(nz)
    if total_vox == 0:
        return {'count': 0, 'placements': [], 'voxel_grid': np.zeros((0,0,0), dtype=bool),
                'grid_origin': bb_min, 'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': None}

    if total_vox > max_voxels:
        raise MemoryError(f"voxel grid too large: shape=({nx},{ny},{nz}), voxels={total_vox}. Reduce resolution or increase max_voxels.")

    # Candidate start index ranges so that cuboid fully lies within bounding box in each axis
    max_start_x = int(np.floor((bb_max[0] - bb_min[0] - Lx) / pitch + 1e-12))
    max_start_y = int(np.floor((bb_max[1] - bb_min[1] - Ly) / pitch + 1e-12))
    max_start_z = int(np.floor((bb_max[2] - bb_min[2] - Lz) / pitch + 1e-12))
    if max_start_x < 0 or max_start_y < 0 or max_start_z < 0:
        # cuboid larger than bounding extents in at least one axis
        return {'count': 0, 'placements': [], 'voxel_grid': np.zeros((nx,ny,nz), dtype=bool),
                'grid_origin': bb_min, 'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': None}

    # Enumerate candidate starts (x0,y0,z0) as integer voxel offsets for cuboid lower corner
    x_starts = np.arange(0, max_start_x + 1, dtype=np.int32)
    y_starts = np.arange(0, max_start_y + 1, dtype=np.int32)
    z_starts = np.arange(0, max_start_z + 1, dtype=np.int32)

    # Build candidate list as flat arrays (we'll filter invalid ones)
    candidate_list = []
    for xi in x_starts:
        for yi in y_starts:
            for zi in z_starts:
                candidate_list.append((int(xi), int(yi), int(zi)))
    if len(candidate_list) == 0:
        return {'count': 0, 'placements': [], 'voxel_grid': np.zeros((nx,ny,nz), dtype=bool),
                'grid_origin': bb_min, 'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': None}

    cand_arr = np.array(candidate_list, dtype=np.int32)  # shape (M,3)
    M = cand_arr.shape[0]

    # For each candidate we will test both orientations (if requested).
    # We'll create a candidate entry per orientation so each candidate-orientation pair is considered.

    # Build lower-corner world coords for candidate starts
    lows_x = bb_min[0] + cand_arr[:,0].astype(float) * pitch
    lows_y = bb_min[1] + cand_arr[:,1].astype(float) * pitch
    lows_z = bb_min[2] + cand_arr[:,2].astype(float) * pitch

    # Build corner points (vectorized) for orientation 0 and optionally orientation 1
    # For a candidate at index m, its 8 corners are:
    # (lx, ly, lz), (lx+Lx, ly, lz), (lx, ly+Ly, lz), (lx, ly, lz+Lz), ... etc

    # Collect corners for all candidate-orientation combinations
    corner_pts = []
    candidate_meta = []  # tuples: (cand_index, orient_id)
    def append_corners_for_orientation(orient_id, Dx, Dy, Dz):
        # corners relative offsets
        offsets = np.array([[0, 0, 0],
                            [Dx, 0, 0],
                            [0, Dy, 0],
                            [0, 0, Dz],
                            [Dx, Dy, 0],
                            [Dx, 0, Dz],
                            [0, Dy, Dz],
                            [Dx, Dy, Dz]], dtype=float)
        # For each candidate, compute 8 corner points
        for m in range(M):
            lx, ly, lz = lows_x[m], lows_y[m], lows_z[m]
            # only add if cuboid stays inside bounding box in this orientation (should be true given start ranges for orientation 0,
            # but for rotated orientation some starts may make cube exceed bounds; we can filter later)
            base = np.array([lx, ly, lz], dtype=float)
            pts8 = base + offsets  # shape (8,3)
            # check bounding quick reject (all coords <= bb_max + tiny)
            # We'll still call mesh.contains for the exact test; here we just store
            corner_pts.append(pts8)
            candidate_meta.append((m, orient_id))
    # orientation 0: (Lx,Ly,Lz)
    append_corners_for_orientation(0, Lx, Ly, Lz)
    # orientation 1: swapped X/Y if requested and distinct
    if allow_z90 and (abs(Lx - Ly) > 1e-12):
        append_corners_for_orientation(1, Ly, Lx, Lz)

    # Flatten corner_pts for a single contains() call
    corner_pts_arr = np.vstack(corner_pts)  # shape (Ncand_orient*8, 3)

    # Use trimesh contains on corners (vectorized). If many points, trimesh is reasonably fast.
    inside_flags = mesh.contains(corner_pts_arr)  # boolean array length N*8

    # Reshape into (N_candidates_orient, 8) and decide validity
    Norient = len(candidate_meta)
    inside_matrix = inside_flags.reshape((Norient, 8))
    valid_mask = inside_matrix.all(axis=1)  # True if all 8 corners inside mesh

    # Build arrays for the valid candidate-orientation pairs
    if not valid_mask.any():
        return {'count': 0, 'placements': [], 'voxel_grid': inside_matrix,  # debugging return
                'grid_origin': bb_min, 'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': None}

    # Build arrays of starts and orientation ids for valid entries
    valid_idx = np.where(valid_mask)[0]
    NV = valid_idx.shape[0]

    starts_x = np.empty(NV, dtype=np.int32)
    starts_y = np.empty(NV, dtype=np.int32)
    starts_z = np.empty(NV, dtype=np.int32)
    orient_ids = np.empty(NV, dtype=np.int32)
    # Precompute footprint area (for ordering heuristic) and voxel index ranges (xi0..xi1) in Python
    footprint = np.empty(NV, dtype=np.int32)
    xi0 = np.empty(NV, dtype=np.int32)
    xi1 = np.empty(NV, dtype=np.int32)
    yi0 = np.empty(NV, dtype=np.int32)
    yi1 = np.empty(NV, dtype=np.int32)
    zi0 = np.empty(NV, dtype=np.int32)
    zi1 = np.empty(NV, dtype=np.int32)

    # Precompute xs,ys,zs arrays as numpy for searchsorted operations
    xs_np = np.array(xs)
    ys_np = np.array(ys)
    zs_np = np.array(zs)

    for idx_out, idx_meta in enumerate(valid_idx):
        cand_index, o_id = candidate_meta[int(idx_meta)]
        sx = cand_arr[cand_index, 0]
        sy = cand_arr[cand_index, 1]
        sz = cand_arr[cand_index, 2]
        starts_x[idx_out] = sx
        starts_y[idx_out] = sy
        starts_z[idx_out] = sz
        orient_ids[idx_out] = o_id
        # dims based on orientation
        if o_id == 0:
            Dx, Dy, Dz = Lx, Ly, Lz
        else:
            Dx, Dy, Dz = Ly, Lx, Lz  # swapped

        # Compute which voxel centers fall within the cuboid [lx, lx+Dx)
        x_low = bb_min[0] + sx * pitch
        x_high = x_low + Dx
        y_low = bb_min[1] + sy * pitch
        y_high = y_low + Dy
        z_low = bb_min[2] + sz * pitch
        z_high = z_low + Dz

        # use searchsorted to get index ranges; centers at bb_min + (i+0.5)*pitch
        xi_start = int(np.searchsorted(xs_np, x_low + 1e-12, side='left'))
        xi_end = int(np.searchsorted(xs_np, x_high - 1e-12, side='right'))  # exclusive
        yi_start = int(np.searchsorted(ys_np, y_low + 1e-12, side='left'))
        yi_end = int(np.searchsorted(ys_np, y_high - 1e-12, side='right'))
        zi_start = int(np.searchsorted(zs_np, z_low + 1e-12, side='left'))
        zi_end = int(np.searchsorted(zs_np, z_high - 1e-12, side='right'))

        # clamp ranges to grid
        xi_start = max(0, min(nx, xi_start))
        xi_end = max(0, min(nx, xi_end))
        yi_start = max(0, min(ny, yi_start))
        yi_end = max(0, min(ny, yi_end))
        zi_start = max(0, min(nz, zi_start))
        zi_end = max(0, min(nz, zi_end))

        xi0[idx_out] = xi_start
        xi1[idx_out] = xi_end
        yi0[idx_out] = yi_start
        yi1[idx_out] = yi_end
        zi0[idx_out] = zi_start
        zi1[idx_out] = zi_end

        # footprint: number of voxel centers in XY plane that the cuboid would cover (speed heuristic)
        footprint[idx_out] = max(1, (xi_end - xi_start) * max(1, (yi_end - yi_start)))

    # Create an ordering of candidates: primary by z (lowest first), then smaller footprint first, then y then x
    # Compute the z coordinate of candidate starts in voxel units (starts_z) already is that
    order_keys = np.empty((NV, 4), dtype=np.int32)
    order_keys[:, 0] = starts_z    # primary: low z first
    order_keys[:, 1] = footprint   # then prefer small footprint
    order_keys[:, 2] = starts_y
    order_keys[:, 3] = starts_x
    # lexsort expects keys in increasing order of significance, so we pass reversed
    order_idx = np.lexsort((order_keys[:, 3], order_keys[:, 2], order_keys[:, 1], order_keys[:, 0]))
    # order_idx is array of indices into 0..NV-1 specifying traversal order

    # For Numba we need occ array shape (nx,ny,nz) boolean
    occ = np.zeros((nx, ny, nz), dtype=np.bool_)

    # Call numba greedy placement: note it expects arrays of length M where M = NV, and order_idx length NV
    placed_local = greedy_place_numba(occ, nx, ny, nz,
                                      xi0, xi1, yi0, yi1, zi0, zi1,
                                      order_idx.astype(np.int32))

    # Build placements list from placed_local flags
    placements = []
    for i in range(NV):
        if placed_local[i]:
            sx = int(starts_x[i]); sy = int(starts_y[i]); sz = int(starts_z[i])
            o_id = int(orient_ids[i])
            if o_id == 0:
                Dx, Dy, Dz = Lx, Ly, Lz
                kernel_vox = (int(np.ceil(Dx / pitch)), int(np.ceil(Dy / pitch)), int(np.ceil(Dz / pitch)))
            else:
                Dx, Dy, Dz = Ly, Lx, Lz
                kernel_vox = (int(np.ceil(Dx / pitch)), int(np.ceil(Dy / pitch)), int(np.ceil(Dz / pitch)))
            # center world coordinate of placed cuboid
            cx = float(bb_min[0] + sx * pitch + Dx / 2.0)
            cy = float(bb_min[1] + sy * pitch + Dy / 2.0)
            cz = float(bb_min[2] + sz * pitch + Dz / 2.0)
            placements.append({'index': (sx, sy, sz),
                               'center_world': (cx, cy, cz),
                               'orientation': o_id,
                               'kernel_voxels': kernel_vox,
                               'dims': (Dx, Dy, Dz)})

    result = {'count': len(placements),
              'placements': placements,
              'voxel_grid': (nx, ny, nz),    # not returning full grid for memory, but grid shape
              'grid_origin': bb_min,
              'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
              'pitch': pitch,
              'mesh': mesh}
    return result

# -------------------------
# Demo main
# -------------------------
if __name__ == '__main__':
    # # build wing
    # wing = make_tapered_wing_box(span=8.0, root_chord=3.0, tip_chord=1.2,
    #                              thickness=0.35, sweep_le=0.4, dihedral_deg=3.0)
    
    # # =============================================================================
    # import sys
    # import pandas as pd
    # wingbox_coordinates = pd.read_csv('wingbox_coordinates.csv')
    # wing = create_mesh_from_coordinates(wingbox_coordinates)
    # # sys.exit('Done.')
    # # =============================================================================
    
    # # =============================================================================
    # from nils.point_loads.system_study.volume.cargocomp_geometry_mesher import build_double_chord_polygon, make_prism_from_polygon
    # # Example parameters
    # r = 2.5
    # theta1_deg = 5.0
    # thetafcs_deg = 45
    # pidiv2minustheta1_deg = 90. - thetafcs_deg   # angle between diameter and first chord
    # pidiv2minustheta2_deg = 90. - theta1_deg   # angle between diameter and second chord
    # theta1 = np.deg2rad(pidiv2minustheta1_deg)
    # theta2 = np.deg2rad(pidiv2minustheta2_deg)
    # height = 15.2
    # n_arc = 80

    # # Build polygon and prism
    # poly = build_double_chord_polygon(r, theta1, theta2, n_arc=n_arc)
    # wing, _, _ = make_prism_from_polygon(poly, height, extrude_direction=+1.0)
    # # =============================================================================
    
    # # =============================================================================
    # from nils.point_loads.system_study.volume.nacelle_geometry_mesher import revolve_profile_to_mesh
    
    # Psi_list_closed = np.array([
    # 0.   , 0.005, 0.01 , 0.015, 0.02 , 0.025, 0.03 , 0.035, 0.04 ,
    # 0.045, 0.05 , 0.055, 0.06 , 0.065, 0.07 , 0.075, 0.08 , 0.085,
    # 0.09 , 0.095, 0.1  , 0.105, 0.11 , 0.115, 0.12 , 0.125, 0.13 ,
    # 0.135, 0.14 , 0.145, 0.15 , 0.155, 0.16 , 0.165, 0.17 , 0.175,
    # 0.18 , 0.185, 0.19 , 0.195, 0.2  , 0.205, 0.21 , 0.215, 0.22 ,
    # 0.225, 0.23 , 0.235, 0.24 , 0.245, 0.25 , 0.255, 0.26 , 0.265,
    # 0.27 , 0.275, 0.28 , 0.285, 0.29 , 0.295, 0.3  , 0.305, 0.31 ,
    # 0.315, 0.32 , 0.325, 0.33 , 0.335, 0.34 , 0.345, 0.35 , 0.355,
    # 0.36 , 0.365, 0.37 , 0.375, 0.38 , 0.385, 0.39 , 0.395, 0.4  ,
    # 0.405, 0.41 , 0.415, 0.42 , 0.425, 0.43 , 0.435, 0.44 , 0.445,
    # 0.45 , 0.455, 0.46 , 0.465, 0.47 , 0.475, 0.48 , 0.485, 0.49 ,
    # 0.495, 0.5  , 0.505, 0.51 , 0.515, 0.52 , 0.525, 0.53 , 0.535,
    # 0.54 , 0.545, 0.55 , 0.555, 0.56 , 0.565, 0.57 , 0.575, 0.58 ,
    # 0.585, 0.59 , 0.595, 0.6  , 0.605, 0.61 , 0.615, 0.62 , 0.625,
    # 0.63 , 0.635, 0.64 , 0.645, 0.65 , 0.655, 0.66 , 0.665, 0.67 ,
    # 0.675, 0.68 , 0.685, 0.69 , 0.695, 0.7  , 0.705, 0.71 , 0.715,
    # 0.72 , 0.725, 0.73 , 0.735, 0.74 , 0.745, 0.75 , 0.755, 0.76 ,
    # 0.765, 0.77 , 0.775, 0.78 , 0.785, 0.79 , 0.795, 0.8  , 0.805,
    # 0.81 , 0.815, 0.82 , 0.825, 0.83 , 0.835, 0.84 , 0.845, 0.85 ,
    # 0.855, 0.86 , 0.865, 0.87 , 0.875, 0.88 , 0.885, 0.89 , 0.895,
    # 0.9  , 0.905, 0.91 , 0.915, 0.92 , 0.925, 0.93 , 0.935, 0.94 ,
    # 0.945, 0.95 , 0.955, 0.96 , 0.965, 0.97 , 0.975, 0.98 , 0.985,
    # 0.99 , 0.995, 1.   ]) * 4

    # zeta_list_closed = np.array([
    # 0.        , 0.02741722, 0.03861335, 0.04705651, 0.05404363,
    # 0.06008184, 0.06543369, 0.07025658, 0.07465406, 0.07869873,
    # 0.08244396, 0.08593042, 0.08919002, 0.09224839, 0.09512649,
    # 0.09784178, 0.100409  , 0.1028407 , 0.10514776, 0.10733962,
    # 0.10942458, 0.11140999, 0.11330238, 0.1151076 , 0.11683093,
    # 0.11847712, 0.12005051, 0.12155505, 0.12299435, 0.12437173,
    # 0.12569026, 0.12695276, 0.12816185, 0.12931996, 0.13042937,
    # 0.13149217, 0.13251034, 0.13348573, 0.13442008, 0.13531502,
    # 0.13617206, 0.13699266, 0.13777818, 0.13852989, 0.13924902,
    # 0.13993671, 0.14059404, 0.14122206, 0.14182173, 0.14239398,
    # 0.1429397 , 0.14345972, 0.14395484, 0.14442581, 0.14487336,
    # 0.14529818, 0.14570091, 0.14608218, 0.14644258, 0.14678269,
    # 0.14710303, 0.14740412, 0.14768645, 0.14795049, 0.14819669,
    # 0.14842546, 0.14863722, 0.14883234, 0.14901119, 0.14917413,
    # 0.14932148, 0.14945357, 0.14957068, 0.14967311, 0.14976113,
    # 0.14983498, 0.14989493, 0.14994119, 0.14997399, 0.14999353,
    # 0.15      , 0.14999312, 0.14997261, 0.14993862, 0.14989133,
    # 0.14983087, 0.14975739, 0.14967102, 0.14957188, 0.14946008,
    # 0.14933573, 0.14919894, 0.14904979, 0.14888836, 0.14871473,
    # 0.14852898, 0.14833115, 0.14812132, 0.14789952, 0.14766581,
    # 0.14742021, 0.14716277, 0.1468935 , 0.14661242, 0.14631955,
    # 0.1460149 , 0.14569846, 0.14537024, 0.14503022, 0.1446784 ,
    # 0.14431476, 0.14393927, 0.14355191, 0.14315264, 0.14274142,
    # 0.14231823, 0.141883  , 0.14143569, 0.14097625, 0.14050461,
    # 0.14002072, 0.13952451, 0.1390159 , 0.13849482, 0.1379612 ,
    # 0.13741495, 0.13685598, 0.13628421, 0.13569953, 0.13510187,
    # 0.1344911 , 0.13386714, 0.13322987, 0.13257918, 0.13191495,
    # 0.13123708, 0.13054544, 0.1298399 , 0.12912034, 0.12838663,
    # 0.12763864, 0.12687623, 0.12609925, 0.12530758, 0.12450106,
    # 0.12367954, 0.12284288, 0.12199092, 0.12112351, 0.12024048,
    # 0.11934168, 0.11842693, 0.11749608, 0.11654896, 0.11558538,
    # 0.11460518, 0.11360818, 0.1125942 , 0.11156306, 0.11051456,
    # 0.10944854, 0.10836478, 0.10726311, 0.10614333, 0.10500524,
    # 0.10384864, 0.10267333, 0.10147911, 0.10026577, 0.0990331 ,
    # 0.0977809 , 0.09650894, 0.09521702, 0.09390492, 0.09257242,
    # 0.09121929, 0.08984532, 0.08845027, 0.08703393, 0.08559606,
    # 0.08413643, 0.0826548 , 0.08115095, 0.07962463, 0.0780756 ,
    # 0.07650362, 0.07490845, 0.07328985, 0.07164756, 0.06998133,
    # 0.06829092, 0.06657607, 0.06483652, 0.06307203, 0.06128233,
    # 0.05946715, 0.05762625, 0.05575935, 0.05386618, 0.05194649,
    # 0.05      ]) * 6

    # # build mesh
    # wing = revolve_profile_to_mesh(Psi_list_closed, zeta_list_closed, n_phi=128//4, cap_ends=True)
    # # =============================================================================
    
    # =============================================================================
    from nils.point_loads.system_study.volume.fcsdisk_geometry_mesher import make_disk
    
    # example parameters
    d = 1.5    # depth in meters
    r = 2.5    # radius in meters

    wing = make_disk(d, r)
    # =============================================================================

    # example 1: medium fine pitch automatic (safe)
    # cuboid = (1.0, 0.5, 0.1)   # (Lx, Ly, Lz) meters
    cuboid = (1.0, 0.5, 0.25)   # (Lx, Ly, Lz) meters
    print("Wing bounds:", wing.bounds)
    print("Cuboid dims (Lx, Ly, Lz):", cuboid)

    # Automatic pitch (safe selection)
    # result = pack_cuboids_in_mesh(wing, cuboid_dims=cuboid, pitch=0.05, allow_z90=True)
    result = pack_cuboids_in_mesh(wing, cuboid_dims=cuboid, pitch=0.1, allow_z90=True)
    print("Pitch used (m):", result['pitch'])
    print("Number of stacks placed:", result['count'])
    for i, p in enumerate(result['placements'][:20]):
        print(f" {i:3d}: center={np.round(p['center_world'],4)}, orient={p['orientation']}, dims={np.round(p['dims'],4)}")

    #%%
    
    # -------------------------
    # Visualization helpers (updated to include rotated dims)
    # -------------------------
    def plot_mesh_and_placements(mesh, pack_result, cuboid_dims, show_voxels=False):
        
        fig = plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111, projection='3d')

        tri_verts = mesh.vertices[mesh.faces]
        mesh_collection = Poly3DCollection(tri_verts, alpha=0.18, linewidths=0.2)
        mesh_collection.set_facecolor((0.4,0.6,0.9,0.18))
        ax.add_collection3d(mesh_collection)

        Lx, Ly, Lz = cuboid_dims
        for p in pack_result['placements']:
            cx, cy, cz = p['center_world']
            orient = p['orientation']
            Dx, Dy, Dz = p['dims']
            hx, hy, hz = Dx/2.0, Dy/2.0, Dz/2.0
            corners = np.array([[cx-hx, cy-hy, cz-hz],
                                [cx+hx, cy-hy, cz-hz],
                                [cx+hx, cy+hy, cz-hz],
                                [cx-hx, cy+hy, cz-hz],
                                [cx-hx, cy-hy, cz+hz],
                                [cx+hx, cy-hy, cz+hz],
                                [cx+hx, cy+hy, cz+hz],
                                [cx-hx, cy+hy, cz+hz]])
            faces_i = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
            face_verts = [[corners[idx] for idx in face] for face in faces_i]
            poly = Poly3DCollection(face_verts, alpha=0.75)
            poly.set_edgecolor((0,0,0,0.6))
            if orient == 0:
                poly.set_facecolor((0.9,0.6,0.3,0.6))
            else:
                poly.set_facecolor((0.3,0.9,0.5,0.6))
            ax.add_collection3d(poly)

        # autoscale
        all_pts = mesh.vertices if not pack_result['placements'] else np.vstack([mesh.vertices] + [[p['center_world']] for p in pack_result['placements']])
        min_v = all_pts.min(axis=0)
        max_v = all_pts.max(axis=0)
        # ax.set_xlim(min_v[0], max_v[0])
        # ax.set_ylim(min_v[1], max_v[1])
        # ax.set_zlim(min_v[2], max_v[2])
        ax.set_aspect('equal')
        ax.set_xlabel('X (chordwise) [m]')
        ax.set_ylabel('Y (spanwise) [m]')
        ax.set_zlabel('Z (vertical) [m]')
        ax.set_title(f"Placed {pack_result['count']} cuboids (pitch={pack_result['pitch']:.4f} m)")
        
        # Remove the grey shading of the coordinate planes and set ticks and ticklabels
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.grid(True)
        ax.xaxis._axinfo['grid']['color'] = (0, 0, 0, 0.1)  # Black grid lines with 0.2 opacity
        ax.yaxis._axinfo['grid']['color'] = (0, 0, 0, 0.1)  # Black grid lines with 0.2 opacity
        ax.zaxis._axinfo['grid']['color'] = (0, 0, 0, 0.1)  # Black grid lines with 0.2 opacity
        # x_ticks = scatter_points[:, 0]
        # y_ticks = scatter_points[:, 1]
        # z_ticks = scatter_points[:, 2]
        # ax.set_xticks(x_ticks)
        # ax.set_yticks(y_ticks)
        # ax.set_zticks(z_ticks)
        # ax.set_xticklabels([f'{val:.2g}' for val in x_ticks])
        # ax.set_yticklabels([f'{val:.2g}' for val in y_ticks])
        # ax.set_zticklabels([f'{val:.2g}' for val in z_ticks])
        
        # ax.view_init(elev=30, azim=225, roll=0)
        ax.view_init(elev=135, azim=-70, roll=-75)
        
        plt.tight_layout()
        plt.show()

    plot_mesh_and_placements(wing, result, cuboid)

    # # Example 2: user-specified pitch (valid case)
    # try:
    #     # This pitch divides Lx exactly (1.0/0.2 = 5) but does NOT divide Ly or Lz.
    #     # The code will raise ValueError because pitch is incompatible with cuboid dims,
    #     # which prevents the partially-outside behaviour you observed.
    #     bad_pitch = 0.2
    #     _ = pack_cuboids_in_mesh(wing, cuboid_dims=cuboid, pitch=bad_pitch, allow_z90=True)
    # except ValueError as e:
    #     print("Correctly rejected bad pitch (example):", e)

    # # Example 3: user provides small-enough pitch explicitly (valid)
    # fine_pitch = 0.02   # small enough relative to min dim 0.1
    # result2 = pack_cuboids_in_mesh(wing, cuboid_dims=cuboid, pitch=fine_pitch, allow_z90=True)
    # print("With fine pitch=0.02, placed:", result2['count'])
