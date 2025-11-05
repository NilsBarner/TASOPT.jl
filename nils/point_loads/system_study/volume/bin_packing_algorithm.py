#!/usr/bin/env python3
"""
wing_box_packer_origin_fixed.py

Fixed packing demo:
 - Avoids use of VoxelGrid.origin (which may not exist).
 - Uses mesh.bounds for voxel-grid mapping.
 - Returns and documents an 'origin' set to the bottom aft edge along the
   fuselage centre line (as requested). This origin is not used for internal
   voxel indexing — it is provided as the user-specified reference point.

Requires:
    pip install numpy scipy trimesh matplotlib

Author: ChatGPT (corrected)
"""
import numpy as np
from scipy.signal import convolve
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# -------------------------
# Geometry helpers
# -------------------------
def make_tapered_wing_box(span=10.0, root_chord=3.0, tip_chord=1.2,
                          thickness=0.3, sweep_le=0.0, dihedral_deg=0.0):
    """
    Build a simple tapered wing-box mesh (two-station prism).
    Coordinate system:
      x = forward (chordwise), y = spanwise, z = up.

    Returns: trimesh.Trimesh (attempts repairs and falls back to convex hull).
    """
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
        # order: top leading, top trailing, bottom trailing, bottom leading
        verts.append([le,       y, z_off + thickness/2.0])
        verts.append([le + c,   y, z_off + thickness/2.0])
        verts.append([le + c,   y, z_off - thickness/2.0])
        verts.append([le,       y, z_off - thickness/2.0])

    verts = np.array(verts)

    def v(i,j): return 4*i + j

    faces = []
    faces += [[v(0,0), v(1,0), v(1,1)], [v(0,0), v(1,1), v(0,1)]]  # top
    faces += [[v(0,3), v(1,3), v(1,2)], [v(0,3), v(1,2), v(0,2)]]  # bottom
    faces += [[v(0,0), v(0,3), v(1,3)], [v(0,0), v(1,3), v(1,0)]]  # leading
    faces += [[v(0,1), v(1,1), v(1,2)], [v(0,1), v(1,2), v(0,2)]]  # trailing

    faces = np.array(faces)
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=True)

    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()   # modifies mesh in-place
        except Exception:
            pass

    if not mesh.is_watertight:
        mesh = trimesh.Trimesh(vertices=verts).convex_hull

    return mesh

# -------------------------
# Packing algorithm
# -------------------------
def pack_cuboids_in_mesh(mesh, cuboid_dims, pitch=None, max_voxels=200**3):
    """
    Pack axis-aligned equal cuboids into the interior of `mesh`.

    Parameters
    ----------
    mesh : trimesh.Trimesh
      Watertight mesh representing the containing volume.
    cuboid_dims : (Lx, Ly, Lz)
      Cuboid dimensions in meters (Lx = chordwise (x), Ly = spanwise (y), Lz = vertical (z)).
    pitch : float or None
      Voxel pitch in meters. If None, chosen from cuboid dims as a heuristic.
    max_voxels : int
      Safety limit for number of voxels.

    Returns
    -------
    result : dict with keys
      - 'count' : int, number of cuboids placed
      - 'placements' : list of dicts {'index':(ix,iy,iz), 'center_world':(x,y,z)}
      - 'voxel_grid' : boolean ndarray (nx,ny,nz) of interior voxels
      - 'grid_origin' : np.array([min_x, min_y, min_z]) used for voxel indexing
      - 'origin' : np.array([aft_x, 0.0, bottom_z])  <-- requested reference point
      - 'pitch' : voxel pitch (m)
      - 'kernel_voxels' : (kx,ky,kz)
      - 'mesh' : original mesh
    Notes
    -----
    - The returned 'origin' is the bottom aft edge along the fuselage centre line
      (computed from mesh bounds) as you requested. It is provided for your reference.
    - Internally, voxel-index -> world-coordinate mapping uses 'grid_origin'
      = mesh.bounds[0] (the minimum corner), i.e. voxel index 0 corresponds to grid_origin.
    """
    Lx, Ly, Lz = cuboid_dims

    # Heuristic pitch selection
    if pitch is None:
        pitch = min(Lx, Ly, Lz) / 3.0

    # integer kernel sizes (>=1)
    kx = max(1, int(round(Lx / pitch)))
    ky = max(1, int(round(Ly / pitch)))
    kz = max(1, int(round(Lz / pitch)))

    # recompute pitch so cuboid maps exactly to integer voxels
    pitch = float(min(Lx / kx if kx else Lx,
                      Ly / ky if ky else Ly,
                      Lz / kz if kz else Lz))

    # voxelize mesh (Trimesh VoxelGrid). We do NOT use vg.origin.
    try:
        vg = mesh.voxelized(pitch)
    except Exception as e:
        raise RuntimeError("Voxelization failed: " + str(e))

    inside = vg.matrix.copy()   # boolean ndarray, shape (nx, ny, nz)
    grid_origin = mesh.bounds[0].astype(float)   # min corner (used for index->world mapping)
    # --- user-requested origin: bottom aft edge on fuselage centre line ---
    # bounds: [[minx, miny, minz], [maxx, maxy, maxz]]
    bb_min, bb_max = mesh.bounds[0], mesh.bounds[1]
    origin_user = np.array([bb_max[0], 0.0, bb_min[2]], dtype=float)
    # -----------------------------------------------------------------------

    nx, ny, nz = inside.shape
    total_vox = int(nx) * int(ny) * int(nz)
    if total_vox > max_voxels:
        raise MemoryError(f"voxel grid too large: shape={inside.shape}, voxels={total_vox}")

    # block-sum via valid convolution
    P = inside.astype(np.int32)
    kernel = np.ones((kx, ky, kz), dtype=np.int32)
    block_sum = convolve(P, kernel, mode='valid')   # shape (nx-kx+1, ny-ky+1, nz-kz+1)
    full = kx * ky * kz
    if block_sum.size == 0:
        return {'count': 0, 'placements': [], 'voxel_grid': inside,
                'grid_origin': grid_origin, 'origin': origin_user,
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': (kx,ky,kz)}

    candidates = (block_sum == full)
    idxs = np.argwhere(candidates)
    if idxs.size == 0:
        return {'count': 0, 'placements': [], 'voxel_grid': inside,
                'grid_origin': grid_origin, 'origin': origin_user,
                'pitch': pitch, 'mesh': mesh, 'kernel_voxels': (kx,ky,kz)}

    # Sort candidates: bottom-first (small z index) then y then x
    order = np.lexsort((idxs[:,0], idxs[:,1], idxs[:,2]))  # gives z primary last -> adjust
    idxs = idxs[order]
    # Ensure z-major ordering (stable)
    idxs = idxs[np.argsort(idxs[:,2], kind='stable')]

    occ = np.zeros_like(inside, dtype=bool)
    placements = []
    for xs, ys, zs in idxs:
        x0, x1 = int(xs), int(xs + kx)
        y0, y1 = int(ys), int(ys + ky)
        z0, z1 = int(zs), int(zs + kz)
        if occ[x0:x1, y0:y1, z0:z1].any():
            continue
        occ[x0:x1, y0:y1, z0:z1] = True
        # Map voxel indices -> world coords using grid_origin (min corner):
        cx = grid_origin[0] + (x0 + kx/2.0) * pitch
        cy = grid_origin[1] + (y0 + ky/2.0) * pitch
        cz = grid_origin[2] + (z0 + kz/2.0) * pitch
        placements.append({'index': (x0, y0, z0), 'center_world': (cx, cy, cz)})

    return {'count': len(placements), 'placements': placements,
            'voxel_grid': inside, 'grid_origin': grid_origin, 'origin': origin_user,
            'pitch': pitch, 'mesh': mesh, 'kernel_voxels': (kx,ky,kz)}

# -------------------------
# Visualization helpers
# -------------------------
def plot_mesh_and_placements(mesh, pack_result, cuboid_dims, show_voxels=False):
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='3d')

    tri_verts = mesh.vertices[mesh.faces]
    mesh_collection = Poly3DCollection(tri_verts, alpha=0.15, linewidths=0.2)
    mesh_collection.set_facecolor((0.5,0.5,0.8,0.15))
    ax.add_collection3d(mesh_collection)

    Lx, Ly, Lz = cuboid_dims
    for p in pack_result['placements']:
        cx, cy, cz = p['center_world']
        hx, hy, hz = Lx/2.0, Ly/2.0, Lz/2.0
        corners = np.array([[cx-hx, cy-hy, cz-hz],
                            [cx+hx, cy-hy, cz-hz],
                            [cx+hx, cy+hy, cz-hz],
                            [cx-hx, cy+hy, cz-hz],
                            [cx-hx, cy-hy, cz+hz],
                            [cx+hx, cy-hy, cz+hz],
                            [cx+hx, cy+hy, cz+hz],
                            [cx-hx, cy+hy, cz+hz] if False else [cx-hx, cy+hy, cz+hz]])  # safe corner expr
        # faces for poly
        faces_i = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
        face_verts = [[corners[idx] for idx in face] for face in faces_i]
        poly = Poly3DCollection(face_verts, alpha=0.6)
        poly.set_edgecolor((0,0,0,0.6))
        poly.set_facecolor((0.9,0.6,0.3,0.3))
        ax.add_collection3d(poly)

    if show_voxels:
        inside = pack_result['voxel_grid']
        grid_origin = pack_result['grid_origin']
        pitch = pack_result['pitch']
        coords = np.argwhere(inside)
        if coords.shape[0] > 10000:
            idx = np.random.choice(coords.shape[0], size=10000, replace=False)
            coords = coords[idx]
        world = grid_origin + (coords + 0.5) * pitch
        ax.scatter(world[:,0], world[:,1], world[:,2], s=1, alpha=0.12)

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
    ax.set_title(f"Placed {pack_result['count']} cuboids")
    plt.tight_layout()
    plt.show()

# -------------------------
# Example main
# -------------------------
if __name__ == '__main__':
    wing = make_tapered_wing_box(span=8.0, root_chord=3.0, tip_chord=1.2,
                                 thickness=0.35, sweep_le=0.4, dihedral_deg=3.0)
    cuboid = (0.1, 0.1, 0.1)

    print("Wing bounds:", wing.bounds)
    print("Cuboid dims (Lx, Ly, Lz):", cuboid)
    result = pack_cuboids_in_mesh(wing, cuboid_dims=cuboid, pitch=None)

    print("Grid origin (min corner for voxel indexing):", result['grid_origin'])
    print("Returned origin (bottom aft edge on centerline):", result['origin'])
    print("Pitch used (m):", result['pitch'])
    print("Kernel voxels (kx,ky,kz):", result.get('kernel_voxels'))
    print("Number of stacks placed:", result['count'])
    for p in result['placements'][:10]:
        print("  ", p['center_world'])

    plot_mesh_and_placements(wing, result, cuboid, show_voxels=False)
