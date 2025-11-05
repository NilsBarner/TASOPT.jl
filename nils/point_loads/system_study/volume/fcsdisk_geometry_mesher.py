#!/usr/bin/env python3
"""
disk_mesh.py

Create a watertight triangular mesh of a 3D disk (a short cylinder) whose
extent in X is `depth` (meters) and whose circular cross-section in the Y-Z
plane has radius `radius` (meters).

Function:
    make_disk(depth, radius) -> trimesh.Trimesh

Only inputs are `depth` and `radius`. The function returns a single, watertight
trimesh.Trimesh object.

Requires:
    pip install numpy trimesh matplotlib   # matplotlib only for the optional demo plot
"""

__all__ = ["make_disk"]

import numpy as np
import trimesh

def make_disk(depth, radius):
    """
    Build and return a watertight mesh representing a disk (short cylinder).

    Parameters
    ----------
    depth : float
        Length in the X direction (m). Must be > 0.
    radius : float
        Radius in the Y-Z plane (m). Must be >= 0.

    Returns
    -------
    mesh : trimesh.Trimesh
        Watertight triangular mesh of the disk (processed).
    """
    if depth <= 0:
        raise ValueError("depth must be > 0")
    if radius < 0:
        raise ValueError("radius must be >= 0")

    # fixed tessellation resolution (good quality & reasonably fast)
    n_phi = 128  # azimuthal subdivisions (must be >= 3)

    half = depth / 2.0
    # =============================================================================
    x_shift = 30
    # =============================================================================
    x0 = -half + x_shift
    x1 = +half + x_shift

    phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    cos = np.cos(phis)
    sin = np.sin(phis)

    verts = []
    # front ring (x = x0)
    for j in range(n_phi):
        verts.append((x0, radius * cos[j], radius * sin[j]))
    # back ring (x = x1)
    for j in range(n_phi):
        verts.append((x1, radius * cos[j], radius * sin[j]))

    # optional centers for caps
    idx_center_front = len(verts); verts.append((x0, 0.0, 0.0))
    idx_center_back  = len(verts); verts.append((x1, 0.0, 0.0))

    faces = []
    # side faces (quads split into two triangles)
    for j in range(n_phi):
        jn = (j + 1) % n_phi
        a = j                # front j
        b = jn               # front j+1
        c = n_phi + jn       # back j+1
        d = n_phi + j        # back j
        faces.append((a, b, c))
        faces.append((a, c, d))

    # front cap (triangulate fan about center_front)
    for j in range(n_phi):
        jn = (j + 1) % n_phi
        # order chosen so normal points towards -X (outward)
        faces.append((idx_center_front, jn, j))

    # back cap (triangulate fan about center_back)
    for j in range(n_phi):
        jn = (j + 1) % n_phi
        # order chosen so normal points towards +X (outward)
        faces.append((idx_center_back, n_phi + j, n_phi + jn))

    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)

    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=True)

    # attempt to fix tiny issues if any (shouldn't be necessary)
    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        # fallback: return convex hull (very unlikely for this simple geometry)
        mesh = trimesh.Trimesh(vertices=verts, faces=faces).convex_hull

    return mesh

# Minimal demo when run directly
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    # example parameters
    d = 1.5    # depth in meters
    r = 2.5    # radius in meters

    m = make_disk(d, r)
    print("Vertices:", m.vertices.shape, "Faces:", m.faces.shape, "Watertight:", m.is_watertight)

    # quick 3D plot for verification
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    tri_verts = m.vertices[m.faces]
    col = Poly3DCollection(tri_verts, alpha=0.9)
    col.set_facecolor((0.6,0.8,1.0,0.9))
    col.set_edgecolor((0.05,0.05,0.05,0.2))
    ax.add_collection3d(col)

    pts = m.vertices
    mins = pts.min(axis=0); maxs = pts.max(axis=0)
    ax.set_xlim(mins[0]-0.05, maxs[0]+0.05)
    ax.set_ylim(mins[1]-r*0.05, maxs[1]+r*0.05)
    ax.set_zlim(mins[2]-r*0.05, maxs[2]+r*0.05)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title(f"Disk (depth={d} m, radius={r} m)")
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()
