#!/usr/bin/env python3
"""
wing_box_packer_numba_prefilter_trimesh_fixed.py

Same packer as you provided, but the voxel prefilter is implemented using
sampling of voxel-centers with trimesh.contains (the approach that matches
the original working routine). This avoids Open3D and fixes the case where
trimesh.voxelized() produced an interior mask incompatible with the
convolution prefilter (which resulted in zero candidates).

This is the minimal-change, self-contained file: only the prefilter function
was changed to use mesh.contains on voxel centers. The rest of your code
(and APIs) were left intact.

Requirements:
    pip install numpy scipy trimesh matplotlib numba

Author: ChatGPT
"""
import time
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
        y = ys[j]; c = chords[j]; le = le_x[j]; z_off = zs[j]
        verts.append([le,       y, z_off + thickness/2.0])
        verts.append([le + c,   y, z_off + thickness/2.0])
        verts.append([le + c,   y, z_off - thickness/2.0])
        verts.append([le,       y, z_off - thickness/2.0])

    def v(i,j): return 4*i + j
    faces = []
    faces += [[v(0,0), v(1,0), v(1,1)], [v(0,0), v(1,1), v(0,1)]]
    faces += [[v(0,3), v(1,3), v(1,2)], [v(0,3), v(1,2), v(0,2)]]
    faces += [[v(0,0), v(0,3), v(1,3)], [v(0,0), v(1,3), v(1,0)]]
    faces += [[v(0,1), v(1,1), v(1,2)], [v(0,1), v(1,2), v(0,2)]]

    mesh = trimesh.Trimesh(vertices=np.array(verts), faces=np.array(faces), process=True)
    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices(); mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        mesh = trimesh.Trimesh(vertices=np.array(verts)).convex_hull
    return mesh

# -------------------------
# Create mesh from uploaded coordinates (kept as in your version)
# -------------------------
import numpy as _np
import trimesh as _trimesh

def create_mesh_from_coordinates(coords):
    xs = coords['xw'].to_numpy()
    ys = coords['yw'].to_numpy()
    zs = coords['zw'].to_numpy()

    verts = _np.vstack((xs, ys, zs)).T

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

    faces = _np.array(faces, dtype=_np.int64)

    mesh = _trimesh.Trimesh(vertices=verts, faces=faces, process=True)

    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        mesh = _trimesh.Trimesh(vertices=verts).convex_hull

    return mesh

# -------------------------
# Numba-accelerated greedy placement (unchanged)
# -------------------------
@njit
def greedy_place_numba(occ, nx, ny, nz,
                       xi0, xi1, yi0, yi1, zi0, zi1,
                       order_idx):
    M = order_idx.shape[0]
    placed = np.zeros(M, dtype=np.bool_)
    for p in range(M):
        cand = order_idx[p]
        x0 = xi0[cand]; x1 = xi1[cand]
        y0 = yi0[cand]; y1 = yi1[cand]
        z0 = zi0[cand]; z1 = zi1[cand]
        if x1 <= x0 or y1 <= y0 or z1 <= z0:
            continue
        overlap = False
        for i in range(x0, x1):
            for j in range(y0, y1):
                for k in range(z0, z1):
                    if occ[i, j, k]:
                        overlap = True
                        break
                if overlap: break
            if overlap: break
        if not overlap:
            placed[cand] = True
            for i in range(x0, x1):
                for j in range(y0, y1):
                    for k in range(z0, z1):
                        occ[i, j, k] = True
    return placed

# -------------------------
# Helper: build voxel grid centers from bounds (unchanged)
# -------------------------
def build_voxel_grid_from_bounds(bb_min, bb_max, pitch):
    xs = np.arange(bb_min[0] + 0.5 * pitch, bb_max[0], pitch)
    ys = np.arange(bb_min[1] + 0.5 * pitch, bb_max[1], pitch)
    zs = np.arange(bb_min[2] + 0.5 * pitch, bb_max[2], pitch)
    return xs, ys, zs

# -------------------------
# Method A: brute corner-check for all candidate starts (unchanged)
# -------------------------
def pack_method_A(mesh, cuboid_dims, pitch, allow_z90):
    Lx, Ly, Lz = cuboid_dims
    bb_min = mesh.bounds[0]; bb_max = mesh.bounds[1]

    max_start_x = int(np.floor((bb_max[0] - bb_min[0] - Lx) / pitch + 1e-12))
    max_start_y = int(np.floor((bb_max[1] - bb_min[1] - Ly) / pitch + 1e-12))
    max_start_z = int(np.floor((bb_max[2] - bb_min[2] - Lz) / pitch + 1e-12))
    if max_start_x < 0 or max_start_y < 0 or max_start_z < 0:
        return {'count':0, 'placements':[], 'time':0.0}

    cand = []
    for xi in range(max_start_x+1):
        for yi in range(max_start_y+1):
            for zi in range(max_start_z+1):
                cand.append((xi, yi, zi))
    M = len(cand)
    corner_pts = []
    meta = []
    for idx, (xi, yi, zi) in enumerate(cand):
        lx = bb_min[0] + xi * pitch
        ly = bb_min[1] + yi * pitch
        lz = bb_min[2] + zi * pitch
        offsets = np.array([[0,0,0],[Lx,0,0],[0,Ly,0],[0,0,Lz],[Lx,Ly,0],[Lx,0,Lz],[0,Ly,Lz],[Lx,Ly,Lz]])
        pts8 = (offsets + np.array([lx,ly,lz])).astype(float)
        corner_pts.append(pts8)
        meta.append((idx,0,xi,yi,zi))
        if allow_z90 and (abs(Lx-Ly) > 1e-12):
            offsets2 = np.array([[0,0,0],[Ly,0,0],[0,Lx,0],[0,0,Lz],[Ly,Lx,0],[Ly,0,Lz],[0,Lx,Lz],[Ly,Lx,Lz]])
            pts8b = (offsets2 + np.array([lx,ly,lz])).astype(float)
            corner_pts.append(pts8b)
            meta.append((idx,1,xi,yi,zi))
    if len(corner_pts) == 0:
        return {'count':0, 'placements':[], 'time':0.0}

    corners_flat = np.vstack(corner_pts)
    t0 = time.perf_counter()
    contains_flags = mesh.contains(corners_flat)
    t1 = time.perf_counter()
    contains_matrix = contains_flags.reshape((len(corner_pts), 8))
    valid_mask = contains_matrix.all(axis=1)
    valid_idx = np.where(valid_mask)[0]
    NV = len(valid_idx)
    if NV == 0:
        return {'count':0, 'placements':[], 'time': t1-t0}

    # prepare arrays for numba greedy
    xs_centers, ys_centers, zs_centers = build_voxel_grid_from_bounds(bb_min, bb_max, pitch)
    xs_np = xs_centers; ys_np = ys_centers; zs_np = zs_centers
    xi0 = np.empty(NV, dtype=np.int32); xi1 = np.empty(NV, dtype=np.int32)
    yi0 = np.empty(NV, dtype=np.int32); yi1 = np.empty(NV, dtype=np.int32)
    zi0 = np.empty(NV, dtype=np.int32); zi1 = np.empty(NV, dtype=np.int32)
    starts_x = np.empty(NV, dtype=np.int32); starts_y = np.empty(NV, dtype=np.int32); starts_z = np.empty(NV, dtype=np.int32)
    orient_ids = np.empty(NV, dtype=np.int32)
    for out_i, k in enumerate(valid_idx):
        idx, orient, xi, yi, zi = meta[k]
        starts_x[out_i]=xi; starts_y[out_i]=yi; starts_z[out_i]=zi; orient_ids[out_i]=orient
        if orient == 0:
            Dx, Dy, Dz = Lx, Ly, Lz
        else:
            Dx, Dy, Dz = Ly, Lx, Lz
        x_low = bb_min[0] + xi * pitch; x_high = x_low + Dx
        y_low = bb_min[1] + yi * pitch; y_high = y_low + Dy
        z_low = bb_min[2] + zi * pitch; z_high = z_low + Dz
        xi_start = int(np.searchsorted(xs_np, x_low + 1e-12, side='left'))
        xi_end   = int(np.searchsorted(xs_np, x_high - 1e-12, side='right'))
        yi_start = int(np.searchsorted(ys_np, y_low + 1e-12, side='left'))
        yi_end   = int(np.searchsorted(ys_np, y_high - 1e-12, side='right'))
        zi_start = int(np.searchsorted(zs_np, z_low + 1e-12, side='left'))
        zi_end   = int(np.searchsorted(zs_np, z_high - 1e-12, side='right'))
        xi0[out_i]=max(0,min(len(xs_np),xi_start)); xi1[out_i]=max(0,min(len(xs_np),xi_end))
        yi0[out_i]=max(0,min(len(ys_np),yi_start)); yi1[out_i]=max(0,min(len(ys_np),yi_end))
        zi0[out_i]=max(0,min(len(zs_np),zi_start)); zi1[out_i]=max(0,min(len(zs_np),zi_end))

    # ordering
    footprint = (xi1-xi0).clip(min=0) * (yi1-yi0).clip(min=0)
    order_keys = np.empty((NV,4), dtype=np.int32)
    order_keys[:,0]=starts_z; order_keys[:,1]=footprint; order_keys[:,2]=starts_y; order_keys[:,3]=starts_x
    order_idx = np.lexsort((order_keys[:,3], order_keys[:,2], order_keys[:,1], order_keys[:,0]))

    occ = np.zeros((len(xs_np), len(ys_np), len(zs_np)), dtype=np.bool_)
    t2 = time.perf_counter()
    placed = greedy_place_numba(occ, occ.shape[0], occ.shape[1], occ.shape[2],
                                xi0, xi1, yi0, yi1, zi0, zi1, order_idx.astype(np.int32))
    t3 = time.perf_counter()
    placements=[]
    for i in range(NV):
        if placed[i]:
            sx = int(starts_x[i]); sy = int(starts_y[i]); sz = int(starts_z[i])
            orient = int(orient_ids[i])
            if orient==0: Dx,Dy,Dz=Lx,Ly,Lz
            else: Dx,Dy,Dz=Ly,Lx,Lz
            cx = float(bb_min[0] + sx*pitch + Dx/2.0)
            cy = float(bb_min[1] + sy*pitch + Dy/2.0)
            cz = float(bb_min[2] + sz*pitch + Dz/2.0)
            placements.append({'index':(sx,sy,sz),'center_world':(cx,cy,cz),'orientation':orient,'dims':(Dx,Dy,Dz)})
    return {'count':len(placements),'placements':placements,'time':(t1-t0)+(t3-t2)}

# -------------------------
# Voxel prefilter using mesh.contains on voxel centers (fixed)
# -------------------------
def voxel_prefilter_trimesh(mesh, pitch):
    """
    Build an inside mask by sampling the voxel centers with mesh.contains.
    Returns:
      - inside (boolean ndarray) of shape (nx,ny,nz) marking interior points
      - xs, ys, zs arrays of voxel-center coordinates
    This mirrors the behavior of the original working code that used mesh.contains
    on grid centers, avoiding the mismatches that can occur with trimesh.voxelized().
    """
    bb_min = mesh.bounds[0].astype(float); bb_max = mesh.bounds[1].astype(float)
    xs = np.arange(bb_min[0] + 0.5*pitch, bb_max[0], pitch)
    ys = np.arange(bb_min[1] + 0.5*pitch, bb_max[1], pitch)
    zs = np.arange(bb_min[2] + 0.5*pitch, bb_max[2], pitch)
    nx, ny, nz = len(xs), len(ys), len(zs)
    if nx == 0 or ny == 0 or nz == 0:
        return np.zeros((nx,ny,nz), dtype=bool), xs, ys, zs
    # sample centers and call mesh.contains once (vectorized)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    pts = np.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
    contains = mesh.contains(pts)
    inside = contains.reshape((nx, ny, nz))
    return inside, xs, ys, zs

# -------------------------
# Method B: voxel prefilter (trimesh contains) -> corner-check -> numba greedy
# -------------------------
def pack_method_B_prefilter(mesh, cuboid_dims, pitch, allow_z90):
    Lx, Ly, Lz = cuboid_dims
    bb_min = mesh.bounds[0]; bb_max = mesh.bounds[1]
    t0 = time.perf_counter()
    acceptable, xs, ys, zs = voxel_prefilter_trimesh(mesh, pitch)
    t1 = time.perf_counter()
    nx, ny, nz = acceptable.shape

    # kernel sizes in voxels for both orientations (rounded)
    kx = max(1, int(round(Lx / pitch))); ky = max(1, int(round(Ly / pitch))); kz = max(1, int(round(Lz / pitch)))

    P = acceptable.astype(np.int32)

    # Full-match first; if none found, relax threshold progressively
    candidates_list = []
    if nx - kx + 1 > 0 and ny - ky + 1 > 0 and nz - kz + 1 > 0:
        block_sum0 = convolve(P, np.ones((kx,ky,kz), dtype=np.int32), mode='valid')
        full_val = kx*ky*kz
        idxs0 = np.argwhere(block_sum0 == full_val)
        if idxs0.size == 0:
            # relax thresholds if no exact match
            thresholds = [0.99, 0.98, 0.95, 0.9, 0.8, 0.6]
            for th in thresholds:
                thr = int(np.floor(full_val * th))
                idxs = np.argwhere(block_sum0 >= thr)
                if idxs.size > 0:
                    idxs0 = idxs
                    break
        for (xs_i, ys_i, zs_i) in idxs0:
            candidates_list.append((int(xs_i), int(ys_i), int(zs_i), 0))

    # rotated orientation
    if allow_z90 and (abs(Lx-Ly) > 1e-12):
        kx2 = ky; ky2 = kx; kz2 = kz
        if nx - kx2 + 1 > 0 and ny - ky2 + 1 > 0 and nz - kz2 + 1 > 0:
            block_sum1 = convolve(P, np.ones((kx2,ky2,kz2), dtype=np.int32), mode='valid')
            full1 = kx2*ky2*kz2
            idxs1 = np.argwhere(block_sum1 == full1)
            if idxs1.size == 0:
                thresholds1 = [0.99, 0.98, 0.95, 0.9, 0.8]
                for th in thresholds1:
                    thr = int(np.floor(full1 * th))
                    idxs_tmp = np.argwhere(block_sum1 >= thr)
                    if idxs_tmp.size > 0:
                        idxs1 = idxs_tmp
                        break
            for (xs_i, ys_i, zs_i) in idxs1:
                candidates_list.append((int(xs_i), int(ys_i), int(zs_i), 1))

    t2 = time.perf_counter()
    if len(candidates_list) == 0:
        return {'count':0, 'placements':[], 'time': (t2-t0)}

    # Build corner points for the reduced candidate set and run mesh.contains once (vectorized)
    corner_pts = []; meta = []
    for (sx, sy, sz, orient) in candidates_list:
        lx = bb_min[0] + sx * pitch; ly = bb_min[1] + sy * pitch; lz = bb_min[2] + sz * pitch
        if orient == 0:
            Dx, Dy, Dz = Lx, Ly, Lz
        else:
            Dx, Dy, Dz = Ly, Lx, Lz
        offsets = np.array([[0,0,0],[Dx,0,0],[0,Dy,0],[0,0,Dz],[Dx,Dy,0],[Dx,0,Dz],[0,Dy,Dz],[Dx,Dy,Dz]])
        pts8 = (offsets + np.array([lx,ly,lz])).astype(float)
        corner_pts.append(pts8)
        meta.append((sx,sy,sz,orient, Dx, Dy, Dz))
    corners_flat = np.vstack(corner_pts)
    t3 = time.perf_counter()
    contains_flags = mesh.contains(corners_flat)
    t4 = time.perf_counter()
    contains_matrix = contains_flags.reshape((len(corner_pts), 8))
    valid_mask = contains_matrix.all(axis=1)
    valid_idx = np.where(valid_mask)[0]
    if valid_idx.size == 0:
        return {'count':0, 'placements':[], 'time': (t4-t0)}
    NV = len(valid_idx)

    # Build arrays for the numba greedy, same as in earlier version
    starts_x = np.empty(NV, dtype=np.int32); starts_y = np.empty(NV, dtype=np.int32); starts_z = np.empty(NV, dtype=np.int32)
    xi0 = np.empty(NV, dtype=np.int32); xi1 = np.empty(NV, dtype=np.int32)
    yi0 = np.empty(NV, dtype=np.int32); yi1 = np.empty(NV, dtype=np.int32)
    zi0 = np.empty(NV, dtype=np.int32); zi1 = np.empty(NV, dtype=np.int32)
    orient_ids = np.empty(NV, dtype=np.int32)
    xs_np = xs; ys_np = ys; zs_np = zs
    for out_i, vi in enumerate(valid_idx):
        sx, sy, sz, orient, Dx, Dy, Dz = meta[vi]
        starts_x[out_i]=sx; starts_y[out_i]=sy; starts_z[out_i]=sz; orient_ids[out_i]=orient
        x_low = bb_min[0] + sx * pitch; x_high = x_low + Dx
        y_low = bb_min[1] + sy * pitch; y_high = y_low + Dy
        z_low = bb_min[2] + sz * pitch; z_high = z_low + Dz
        xi_start = int(np.searchsorted(xs_np, x_low + 1e-12, side='left'))
        xi_end   = int(np.searchsorted(xs_np, x_high - 1e-12, side='right'))
        yi_start = int(np.searchsorted(ys_np, y_low + 1e-12, side='left'))
        yi_end   = int(np.searchsorted(ys_np, y_high - 1e-12, side='right'))
        zi_start = int(np.searchsorted(zs_np, z_low + 1e-12, side='left'))
        zi_end   = int(np.searchsorted(zs_np, z_high - 1e-12, side='right'))
        xi0[out_i]=max(0,min(nx,xi_start)); xi1[out_i]=max(0,min(nx,xi_end))
        yi0[out_i]=max(0,min(ny,yi_start)); yi1[out_i]=max(0,min(ny,yi_end))
        zi0[out_i]=max(0,min(nz,zi_start)); zi1[out_i]=max(0,min(nz,zi_end))

    footprint = (xi1-xi0).clip(min=0)*(yi1-yi0).clip(min=0)
    order_keys = np.empty((NV,4), dtype=np.int32)
    order_keys[:,0]=starts_z; order_keys[:,1]=footprint; order_keys[:,2]=starts_y; order_keys[:,3]=starts_x
    order_idx = np.lexsort((order_keys[:,3], order_keys[:,2], order_keys[:,1], order_keys[:,0]))

    occ = np.zeros((nx, ny, nz), dtype=np.bool_)
    t5 = time.perf_counter()
    placed = greedy_place_numba(occ, nx, ny, nz,
                                xi0, xi1, yi0, yi1, zi0, zi1,
                                order_idx.astype(np.int32))
    t6 = time.perf_counter()
    placements=[]
    for i in range(NV):
        if placed[i]:
            sx = int(starts_x[i]); sy = int(starts_y[i]); sz = int(starts_z[i]); orient=int(orient_ids[i])
            Dx = meta[valid_idx[i]][4]; Dy = meta[valid_idx[i]][5]; Dz = meta[valid_idx[i]][6]
            cx = float(bb_min[0] + sx*pitch + Dx/2.0)
            cy = float(bb_min[1] + sy*pitch + Dy/2.0)
            cz = float(bb_min[2] + sz*pitch + Dz/2.0)
            placements.append({'index':(sx,sy,sz),'center_world':(cx,cy,cz),'orientation':orient,'dims':(Dx,Dy,Dz)})
    return {'count':len(placements),'placements':placements,'time':(t1-t0)+(t2-t1)+(t4-t3)+(t6-t5)}

# -------------------------
# Simple plotting helper (unchanged)
# -------------------------
def plot_mesh_and_placements(mesh, pack_result, cuboid_dims):
    fig = plt.figure(figsize=(10,7)); ax = fig.add_subplot(111, projection='3d')
    tri_verts = mesh.vertices[mesh.faces]
    mesh_col = Poly3DCollection(tri_verts, alpha=0.12); mesh_col.set_facecolor((0.4,0.6,0.9,0.12)); ax.add_collection3d(mesh_col)
    for p in pack_result['placements']:
        cx,cy,cz = p['center_world']; Dx,Dy,Dz = p['dims']
        corners = np.array([[cx-Dx/2,cy-Dy/2,cz-Dz/2],
                            [cx+Dx/2,cy-Dy/2,cz-Dz/2],
                            [cx+Dx/2,cy+Dy/2,cz-Dz/2],
                            [cx-Dx/2,cy+Dy/2,cz-Dz/2],
                            [cx-Dx/2,cy-Dy/2,cz+Dz/2],
                            [cx+Dx/2,cy-Dy/2,cz+Dz/2],
                            [cx+Dx/2,cy+Dy/2,cz+Dz/2],
                            [cx-Dx/2,cy+Dy/2,cz+Dz/2]])
        faces_i = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
        face_verts = [[corners[idx] for idx in face] for face in faces_i]
        poly = Poly3DCollection(face_verts, alpha=0.7)
        poly.set_facecolor((0.9,0.6,0.3,0.7)); poly.set_edgecolor((0,0,0,0.5)); ax.add_collection3d(poly)
    mins = mesh.vertices.min(axis=0); maxs = mesh.vertices.max(axis=0)
    ax.set_xlim(mins[0], maxs[0]); ax.set_ylim(mins[1], maxs[1]); ax.set_zlim(mins[2], maxs[2])
    plt.title(f"Placed {pack_result['count']} cuboids"); plt.show()

# -------------------------
# Demo & benchmark (uses the trimesh prefilter implementation)
# -------------------------
if __name__ == '__main__':
    # =============================================================================
    import sys
    import pandas as pd
    wingbox_coordinates = pd.read_csv('wingbox_coordinates.csv')
    wing = create_mesh_from_coordinates(wingbox_coordinates)
    # sys.exit('Done.')
    # =============================================================================

    # example 1: medium fine pitch automatic (safe)
    cuboid = (1.0, 0.5, 0.5)   # (Lx, Ly, Lz) meters

    # auto_pitch = min(cuboid) / 3.0
    
    tA0 = time.perf_counter()
    resA = pack_method_A(wing, cuboid, pitch=0.05, allow_z90=True)
    tA = time.perf_counter()-tA0
    print(f"Method A: placed {resA['count']} in {resA['time']:.3f}s (total elapsed {tA:.3f}s)")

    # Run Method B but using the trimesh voxel prefilter (mesh.contains on centers)
    tB0 = time.perf_counter()
    resB = pack_method_B_prefilter(wing, cuboid, pitch=0.05, allow_z90=True)
    tB = time.perf_counter()-tB0
    print(f"Method B (trimesh voxel prefilter): placed {resB['count']} in (prefilter+checks+greedy) {resB['time']:.3f}s (total elapsed {tB:.3f}s)")
    plot_mesh_and_placements(wing, resB, cuboid)
