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
import sys
import pandas as pd
import numpy as np
import trimesh
import matplotlib.pyplot as plt
from scipy.signal import convolve
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numba import njit

from nils.point_loads.system_study.volume.wingbox_geometry_mesher import create_mesh_from_coordinates
from nils.point_loads.system_study.volume.cargocomp_geometry_mesher import (
    build_double_chord_polygon,
    make_prism_from_polygon,
)
from streamlined_body_generator import generate_streamlined_body_geometry, find_beta_threshold
from nils.point_loads.system_study.volume.nacelle_geometry_mesher import revolve_profile_to_mesh
from nils.point_loads.system_study.volume.fcsdisk_geometry_mesher import make_disk

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
def pack_cuboids_in_mesh(
    mesh, cuboid_dims, pitch=None,
    allow_z90=False, allow_x90=False, allow_y90=False,
    max_voxels=200**3, pitch_tol=1e-8
):

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

    # # Candidate start index ranges so that cuboid fully lies within bounding box in each axis
    # # ensure candidate range is large enough for any rotated orientation:
    # max_dim = max(Lx, Ly, Lz)
    # max_start_x = int(np.floor((bb_max[0] - bb_min[0] - max_dim) / pitch + 1e-12))
    # max_start_y = int(np.floor((bb_max[1] - bb_min[1] - max_dim) / pitch + 1e-12))
    # max_start_z = int(np.floor((bb_max[2] - bb_min[2] - max_dim) / pitch + 1e-12))
    # if max_start_x < 0 or max_start_y < 0 or max_start_z < 0:
    #     # cuboid larger than bounding extents in at least one axis
    #     return {'count': 0, 'placements': [], 'voxel_grid': np.zeros((nx,ny,nz), dtype=bool),
    #             'grid_origin': bb_min, 'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
    #             'pitch': pitch, 'mesh': mesh, 'kernel_voxels': None}
    
    # =============================================================================
    # Candidate start index ranges so that cuboid fully lies within bounding box in each axis
    max_start_x = int(np.floor((bb_max[0] - bb_min[0] - Lx) / pitch + 1e-12))
    max_start_y = int(np.floor((bb_max[1] - bb_min[1] - Ly) / pitch + 1e-12))
    max_start_z = int(np.floor((bb_max[2] - bb_min[2] - Lz) / pitch + 1e-12))
    if max_start_x < 0 or max_start_y < 0 or max_start_z < 0:
        # cuboid larger than bounding extents in at least one axis
        return {'count': 0, 'placements': [], 'voxel_grid': np.zeros((nx,ny,nz), dtype=bool),
                'grid_origin': bb_min, 'origin': np.array([bb_max[0], 0.0, bb_min[2]]),
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': None}
    # =============================================================================

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
    
    # --- build orientation list (closure under requested 90° swaps) ---
    # base orientation
    orient_list = [(float(Lx), float(Ly), float(Lz))]
    
    # allowed transpositions corresponding to 90deg rotations about axes:
    swaps = []
    if allow_z90:
        swaps.append((0, 1))   # swap X <-> Y (rotation about Z)
    if allow_x90:
        swaps.append((1, 2))   # swap Y <-> Z (rotation about X)
    if allow_y90:
        swaps.append((0, 2))   # swap X <-> Z (rotation about Y)
    
    # closure: repeatedly apply allowed swaps to generate all reachable permutations
    i = 0
    while i < len(orient_list):
        cur = orient_list[i]
        for (a, b) in swaps:
            tmp = list(cur)
            tmp[a], tmp[b] = tmp[b], tmp[a]
            tup = (float(tmp[0]), float(tmp[1]), float(tmp[2]))
            if tup not in orient_list:
                orient_list.append(tup)
        i += 1

    # now append corners for every orientation in orient_list
    for oid, (Dx, Dy, Dz) in enumerate(orient_list):
        append_corners_for_orientation(oid, Dx, Dy, Dz)

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
        # dims based on orientation (lookup from orient_list)
        Dx, Dy, Dz = orient_list[int(o_id)]

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
            Dx, Dy, Dz = orient_list[o_id]
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
    
    mesh_list = []
    
    # Wing box
    wingbox_coordinates = pd.read_csv('wingbox_coordinates.csv')
    mesh = create_mesh_from_coordinates(wingbox_coordinates)
    mesh_list.append(mesh)
    
    # Under cabin floor
    
    r = 2.5
    theta1_deg = 5.0
    thetafcs_deg = 30
    pidiv2minustheta1_deg = 90. - thetafcs_deg   # angle between diameter and first chord
    pidiv2minustheta2_deg = 90. - theta1_deg   # angle between diameter and second chord
    theta1 = np.deg2rad(pidiv2minustheta1_deg)
    theta2 = np.deg2rad(pidiv2minustheta2_deg)
    height = 15.2
    n_arc = 80

    poly = build_double_chord_polygon(r, theta1, theta2, n_arc=n_arc)
    mesh, _, _ = make_prism_from_polygon(poly, height, extrude_direction=+1.0)
    mesh_list.append(mesh)
    
    # Nacelle
    
    R_le_over_c = 0.06
    Psi_zeta_max = 0.4
    zeta_max = 2.0/5.0/2.0   # 0.1
    zeta_te = 0.0
    last_ok, first_bad = find_beta_threshold(R_le_over_c, Psi_zeta_max, zeta_max, zeta_te, step_deg=0.01)

    out_dict = generate_streamlined_body_geometry(
        R_le_over_c=zeta_max / 2,
        beta_tail=last_ok/2,
        Psi_zeta_max=Psi_zeta_max,
        zeta_max=zeta_max,
        zeta_te=0.05,
        dimensional_known_dict={'length': 4}
    )
    Psi_list_closed = out_dict['Psi'] * 5
    zeta_list_closed = out_dict['zeta'] * 5
    
    mesh = revolve_profile_to_mesh(Psi_list_closed, zeta_list_closed, n_phi=128//4, cap_ends=True)
    mesh_list.append(mesh)
    
    # Behind cabin
    d = 1.5    # depth in meters
    r = 2.5    # radius in meters
    mesh = make_disk(d, r)
    mesh_list.append(mesh)

    # example 1: medium fine pitch automatic (safe)
    # cuboid = (1.0, 0.5, 0.1)   # (Lx, Ly, Lz) meters
    cuboid = (1.0, 0.5, 0.25)   # (Lx, Ly, Lz) meters
    print("Wing bounds:", mesh.bounds)
    print("Cuboid dims (Lx, Ly, Lz):", cuboid)
    
    # Automatic pitch (safe selection)
    # result = pack_cuboids_in_mesh(mesh, cuboid_dims=cuboid, pitch=0.05, allow_z90=True)
    # result = pack_cuboids_in_mesh(mesh, cuboid_dims=cuboid, pitch=0.1, allow_z90=True)
    
    result_keys = ["wing", "underfloor", "nacelle", "aft"]
    
    result_values = [0, 0, 0, 0]
    volume_ratios = []
    for i, mesh in enumerate(mesh_list):
        
        # if i == 1:
        #     continue
    
        result = pack_cuboids_in_mesh(mesh, cuboid_dims=cuboid, pitch=0.1, allow_z90=False, allow_x90=False, allow_y90=False)
        print("Pitch used (m):", result['pitch'])
        print("Number of stacks placed:", result['count'])
        result_values[i] = result
        
        # <<< ADDED: compute volumes and ratio, print and store
        # single cuboid volume
        single_cuboid_vol = float(np.prod(cuboid))
        # total cuboids volume placed
        total_cuboids_vol = result['count'] * single_cuboid_vol
        # mesh (fluid) volume - use absolute value in case of signed volume
        mesh_vol = abs(float(getattr(mesh, "volume", 0.0)))
        if mesh_vol <= 0.0:
            print(f"Warning: mesh index {i} has non-positive volume ({mesh_vol}). Ratio set to NaN.")
            ratio = np.nan
        else:
            ratio = total_cuboids_vol / mesh_vol
        print(f" Mesh #{i}: mesh_volume = {mesh_vol:.6f} m^3, cuboids_total_volume = {total_cuboids_vol:.6f} m^3, ratio = {ratio:.6f}")
        volume_ratios.append(ratio)
        # <<< /ADDED

    result_dict = dict(zip(result_keys, result_values))
    
    
    # sys.exit('Done.')

    #%%
    
    import pyvista as pv
    
    def trimesh_to_pv(mesh: trimesh.Trimesh) -> pv.PolyData:
        """Convert a trimesh.Trimesh (triangulated) to pyvista.PolyData robustly."""
        verts = mesh.vertices.astype(np.float64)
        faces = mesh.faces.astype(np.int64)
        # VTK faces format: [nVerts, v0, v1, v2, nVerts, v0, v1, v2, ...]
        n_faces = faces.shape[0]
        faces_vtk = np.hstack([np.full((n_faces, 1), 3, dtype=np.int64), faces]).astype(np.int64)
        faces_flat = faces_vtk.ravel()
        return pv.PolyData(verts, faces_flat)
    
    import copy
    from matplotlib_custom_settings import *
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    
    # -------------------------
    # Visualization helpers (updated to include rotated dims)
    # -------------------------
    def plot_mesh_and_placements(
        mesh_list, result_values, cuboid_dims, show_voxels=False
    ):
        
        fig = plt.figure(figsize=(11,8))
        ax = fig.add_subplot(111, projection='3d')
        
        _colors = [colors[1], colors[0], colors[2], colors[3]]

        for i, (mesh, pack_result) in enumerate(zip(mesh_list, result_values)):

            tri_verts = mesh.vertices[mesh.faces]
            mesh_collection = Poly3DCollection(tri_verts, alpha=0.25, linewidths=0.2)
            mesh_collection.set_facecolor(_colors[i])
            # mesh_collection.set_edgecolor("black")
            ax.add_collection3d(mesh_collection)
        
            # if i == 1:
            #     continue
            
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
                poly.set_edgecolor((0,0,0,0.5))
                # if orient == 0:
                #     poly.set_facecolor((0.9,0.6,0.3,0.6))
                # else:
                #     poly.set_facecolor((0.3,0.9,0.5,0.6))
                poly.set_facecolor(_colors[i])
                poly.set_alpha(0.5)
                ax.add_collection3d(poly)

        mesh_mirrored_1 = copy.deepcopy(mesh_list[0])
        mesh_mirrored_1.vertices[:, 1] *= -1
        tri_verts = mesh_mirrored_1.vertices[mesh_mirrored_1.faces]
        mesh_collection = Poly3DCollection(tri_verts, alpha=0.25, linewidths=0.2)
        mesh_collection.set_facecolor(colors[1])
        ax.add_collection3d(mesh_collection)
        
        mesh_mirrored_2 = copy.deepcopy(mesh_list[2])
        mesh_mirrored_2.vertices[:, 1] *= -1
        tri_verts = mesh_mirrored_2.vertices[mesh_mirrored_2.faces]
        mesh_collection = Poly3DCollection(tri_verts, alpha=0.25, linewidths=0.2)
        mesh_collection.set_facecolor(colors[2])
        ax.add_collection3d(mesh_collection)
        
        __colors = [colors[1], colors[2]]
        for i, pack_result in enumerate([result_values[0], result_values[2]]):
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
                corners[:,1] *= -1
                faces_i = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
                face_verts = [[corners[idx] for idx in face] for face in faces_i]
                poly = Poly3DCollection(face_verts, alpha=0.75)
                poly.set_edgecolor((0,0,0,0.5))
                # if orient == 0:
                #     poly.set_facecolor((0.9,0.6,0.3,0.6))
                # else:
                #     poly.set_facecolor((0.3,0.9,0.5,0.6))
                poly.set_facecolor(__colors[i])
                poly.set_alpha(0.5)
                ax.add_collection3d(poly)

        # # autoscale
        # all_pts = mesh.vertices if not pack_result['placements'] else np.vstack([mesh.vertices] + [[p['center_world']] for p in pack_result['placements']])
        # min_v = all_pts.min(axis=0)
        # max_v = all_pts.max(axis=0)
        # ax.set_xlim(min_v[0], max_v[0])
        # ax.set_ylim(min_v[1], max_v[1])
        # ax.set_zlim(min_v[2], max_v[2])
        
        # =============================================================================
        print(np.hstack((np.array(mesh_list), mesh_mirrored_1, mesh_mirrored_2)).flatten())
        for mesh in np.hstack((np.array(mesh_list), mesh_mirrored_1, mesh_mirrored_2)).flatten():
            pv_mesh = pv.wrap(trimesh_to_pv(mesh))  # convert trimesh -> pyvista (see earlier snippet)
            # Extract feature edges (sharp edges). feature_angle is in degrees.
            edges = pv_mesh.extract_feature_edges(feature_angle=30.0,
                                                  boundary_edges=False,
                                                  feature_edges=True,
                                                  non_manifold_edges=False)
            
            # Get the line connectivity (pairs of vertex indices)
            lines = edges.lines.reshape(-1, 3)[:, 1:]  # shape (n_lines, 2)
            points = edges.points                      # shape (n_points, 3)
            
            # Convert into line segments for Matplotlib
            segments = np.array([[points[i], points[j]] for i, j in lines])
            
            lc = Line3DCollection(segments, colors='k', linewidths=1.0, alpha=0.5)
            ax.add_collection3d(lc)
        # =============================================================================
        
        # # =============================================================================
        # # Collect all vertex coordinates from every mesh
        # all_vertices = np.vstack([mesh.vertices for mesh in mesh_list])
        
        # # Compute the overall bounding box
        # mins = all_vertices.min(axis=0)
        # maxs = all_vertices.max(axis=0)
        
        # # Add a small margin (optional)
        # margin = 0.05 * np.linalg.norm(maxs - mins)
        # ax.set_xlim(mins[0] - margin, maxs[0] + margin)
        # ax.set_ylim(mins[1] - margin, maxs[1] + margin)
        # ax.set_zlim(mins[2] - margin, maxs[2] + margin)
        
        # # Force equal aspect
        # ax.set_box_aspect((maxs - mins))  # same scale for all axes
        
        # # ax.set_axis_off()
        # plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        # # =============================================================================
        
        ax.set_aspect('equal')
        ax.set_xlabel('x (m)', labelpad=20)
        ax.set_ylabel('y (m)', labelpad=40)
        ax.set_zlabel('z (m)', labelpad=30)
        # ax.set_title(f"Placed {pack_result['count']} cuboids (pitch={pack_result['pitch']:.4f} m)")
        
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
        # ax.view_init(elev=135, azim=-70, roll=-75)
        ax.view_init(elev=30, azim=225, roll=0)
        
        ax.set_axis_off()
        
        plt.tight_layout()
        plt.show()

    plot_mesh_and_placements(mesh_list, result_values, cuboid)
    
    
    
    sys.exit("What comes after this is not being used.")
    
    #%%
    
    #!/usr/bin/env python3
    """
    Minimal PyVista equivalent of Matplotlib 3D mesh plotting used earlier.
    
    - Converts trimesh.Trimesh -> pyvista.PolyData (robust).
    - Adds transparent faces + black edges.
    - Adds cuboid placements (pv.Cube) from packer results.
    - Produces a tight, cropped interactive scene with no extra white border.
    
    Usage:
        Provide `mesh_list` (list of trimesh.Trimesh) and `result_values`
        (list of pack-result dicts). If not present, demo geometry is used.
    """
    import numpy as np
    import trimesh
    import pyvista as pv
    
    
    def trimesh_to_pv(mesh: trimesh.Trimesh) -> pv.PolyData:
        """Convert a trimesh.Trimesh (triangulated) to pyvista.PolyData robustly."""
        verts = mesh.vertices.astype(np.float64)
        faces = mesh.faces.astype(np.int64)
        # VTK faces format: [nVerts, v0, v1, v2, nVerts, v0, v1, v2, ...]
        n_faces = faces.shape[0]
        faces_vtk = np.hstack([np.full((n_faces, 1), 3, dtype=np.int64), faces]).astype(np.int64)
        faces_flat = faces_vtk.ravel()
        return pv.PolyData(verts, faces_flat)
    
    
    def plot_with_pyvista(mesh_list, result_values=None, colors=None, show=True):
        """
        Plot meshes and placements using PyVista.
        - mesh_list: list of trimesh.Trimesh or pyvista.PolyData
        - result_values: list of dicts; each dict should have key 'placements' (list of placement dicts)
          with fields 'center_world' (tuple), 'dims'=(Dx,Dy,Dz), and 'orientation' (optional)
        - colors: list of RGBA tuples or matplotlib-like color strings
        """
        # Convert inputs
        pv_meshes = []
        for m in mesh_list:
            if isinstance(m, pv.PolyData):
                pv_meshes.append(m)
            elif isinstance(m, trimesh.Trimesh):
                pv_meshes.append(trimesh_to_pv(m))
            else:
                raise TypeError("mesh_list elements must be trimesh.Trimesh or pyvista.PolyData")
    
        if colors is None:
            colors = [(0.4, 0.6, 0.9, 0.18), (0.9, 0.6, 0.3, 0.18), (0.3, 0.9, 0.5, 0.18), (0.9, 0.3, 0.6, 0.18)]
    
        p = pv.Plotter(window_size=(1200, 800), notebook=False, off_screen=not show)
    
        # Add each mesh: faces transparent, black edges (show_edges True)
        for i, pm in enumerate(pv_meshes):
            c = colors[i % len(colors)]
            # pyvista expects RGB (no alpha) for color; control opacity via 'opacity'
            rgb = c[:3] if len(c) >= 3 else c
            alpha = c[3] if len(c) >= 4 else 1.0
            p.add_mesh(pm, color=rgb, opacity=float(alpha), show_edges=True, edge_color='black', lighting=True)
    
        # Add cuboid placements (if provided)
        if result_values is not None:
            for rv_idx, rv in enumerate(result_values):
                placements = rv.get('placements', [])
                for pl in placements:
                    cx, cy, cz = pl['center_world']
                    Dx, Dy, Dz = pl['dims']
                    # create cube: PyVista's Cube uses center and x_length etc
                    cube = pv.Cube(center=(cx, cy, cz), x_length=Dx, y_length=Dy, z_length=Dz)
                    # color cubes by orientation if available
                    orient = pl.get('orientation', 0)
                    cube_color = (0.9, 0.6, 0.3) if orient == 0 else (0.3, 0.9, 0.5)
                    p.add_mesh(cube, color=cube_color, opacity=0.75, show_edges=True, edge_color='k')
    
        # Tight framing: compute combined bounds and set camera / tight view
        # Combine bounds of all pv_meshes and placement cubes
        all_bounds = np.array([pm.bounds for pm in pv_meshes])  # shape (N, 6): (xmin,xmax,ymin,ymax,zmin,zmax)
        xmin = np.min(all_bounds[:, 0])
        xmax = np.max(all_bounds[:, 1])
        ymin = np.min(all_bounds[:, 2])
        ymax = np.max(all_bounds[:, 3])
        zmin = np.min(all_bounds[:, 4])
        zmax = np.max(all_bounds[:, 5])
    
        # Expand slightly for margin
        diag = np.linalg.norm([xmax - xmin, ymax - ymin, zmax - zmin])
        margin = 0.02 * max(diag, 1.0)
        xmin -= margin; xmax += margin
        ymin -= margin; ymax += margin
        zmin -= margin; zmax += margin
    
        # # # Set camera to nicely frame the scene
        # # p.set_background('white')
        # # p.remove_scalar_bar()
        # # p.show_bounds(grid='front', all_edges=True, location='outer', color='black', ticks=(5,5,5))
    
        # # Use the computed bounds to reset camera and tight-fit
        # p.camera_position = 'iso'
        # # Force the camera to look at the center of bounds and set a distance proportional to diag
        # center = ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0)
        # p.camera_set = True
        # p.camera.focal_point = center
        # # Place camera at a sensible offset along current cam direction
        # cam_dir = np.array(p.camera.position) - np.array(p.camera.focal_point)
        # cam_dir = cam_dir / (np.linalg.norm(cam_dir) + 1e-12)
        # p.camera.position = tuple(np.array(center) + cam_dir * (diag * 1.2 + 1e-6))
    
        # Disable axes/labels (you can enable if you like)
        p.show_axes()          # optional: comment out to hide axes
        # p.hide_axes()        # uncomment to hide axes completely
    
        # Lighting and style
        p.enable_eye_dome_lighting()      # improves perception for transparent surfaces
    
        if show:
            p.show()
        else:
            # Offscreen render to an image (example)
            img = p.screenshot(return_img=True)
            return img
    
    
    # if __name__ == "__main__":
    # Demo: if you don't have mesh_list/result_values provide, make simple demo geometry
    demo_meshes = []
    # long box (like your wingbox)
    box = trimesh.creation.box(extents=(12.0, 4.0, 1.0))
    # translate to a position similar to your earlier examples
    box.apply_translation((18.0, 2.0, 0.0))
    demo_meshes.append(box)
    # cone/skewed plate (use triangular prism for demo)
    prism = trimesh.creation.box(extents=(1.0, 3.5, 8.0))
    prism.apply_translation((20.0, 7.0, 1.0))
    demo_meshes.append(prism)
    # cylinder cap at the trailing edge
    cyl = trimesh.creation.cylinder(radius=1.9, height=1.2, sections=64)
    cyl.apply_translation((30.5, 8.0, 0.0))
    demo_meshes.append(cyl)

    # Demo placements: a few cuboids
    demo_results = [{'placements': [
        {'center_world': (20.3554, 2.1034, 0.3), 'dims': (1.0, 0.5, 0.5), 'orientation': 0},
        {'center_world': (21.6554, 2.6034, 0.3), 'dims': (1.0, 0.5, 0.5), 'orientation': 1},
        {'center_world': (24.3554, 2.1034, 0.3), 'dims': (1.0, 0.5, 0.5), 'orientation': 0},
    ]}]

    plot_with_pyvista(demo_meshes, result_values=demo_results)


    