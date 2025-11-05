__all__ = ["revolve_profile_to_mesh"]

import numpy as np
import matplotlib.pyplot as plt
"""
from streamlined_body_generator import generate_streamlined_body_geometry, find_beta_threshold

R_le_over_c = 0.06
Psi_zeta_max = 0.4
zeta_max = 1.5/5.0/2.0   # 0.1
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
Psi_list_closed = out_dict['Psi']
zeta_list_closed = out_dict['zeta']

# Plot below-threshold geometry (your plotting snippet)
Psi = out_dict['Psi']
zeta = out_dict['zeta']
fig, ax = plt.subplots()
ax.plot(Psi, zeta, color='blue', label=f'beta={last_ok:.2f} (below thresh)', linestyle='dashed')
ax.set_xlabel('Chord (m)')
ax.set_ylabel('Thickness (m)')
ax.set_aspect('equal')
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='x', which='both', top=False)
ax.tick_params(axis='y', which='both', right=False)
ax.tick_params(bottom=False)
plt.show()
"""
#%%

#!/usr/bin/env python3
"""
revolution_from_profile.py

Create a watertight triangular mesh by revolving a 2D profile (x = Psi, r = zeta)
around the X-axis (body of revolution). Returns a trimesh.Trimesh object.

Function:
    revolve_profile_to_mesh(x_profile, r_profile, n_phi=96, cap_ends=True, tol_zero=1e-12)

If run as script, the provided Psi_list_closed and zeta_list_closed arrays are used
to build and plot the mesh for visual verification.

Requires:
    pip install numpy trimesh matplotlib

Author: ChatGPT
"""

import sys
from typing import Tuple
import numpy as np
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def revolve_profile_to_mesh(x_profile: np.ndarray,
                            r_profile: np.ndarray,
                            n_phi: int = 96,
                            cap_ends: bool = True,
                            tol_zero: float = 1e-12,
                            process: bool = True) -> trimesh.Trimesh:
    """
    Revolve a profile (x_profile, r_profile) around the X-axis to make a watertight mesh.

    Inputs:
      x_profile : 1D array of x coordinates (length N)
      r_profile : 1D array of radii (>=0) at those x coordinates (length N)
      n_phi     : number of azimuthal samples (>=3)
      cap_ends  : whether to create circular caps at the two ends when radius>tol_zero
      tol_zero  : treated as zero radius if r <= tol_zero (used to create center vertex)
      process   : pass to trimesh.Trimesh(..., process=process)

    Returns:
      trimesh.Trimesh (watertight if caps True and profile sensible)
    """
    x = np.asarray(x_profile, dtype=float).ravel()
    r = np.asarray(r_profile, dtype=float).ravel()
    if x.shape != r.shape:
        raise ValueError("x_profile and r_profile must have same shape")
    if x.size < 2:
        raise ValueError("profile must contain at least two points")

    N = x.size
    M = int(max(3, int(n_phi)))

    # precompute azimuthal angles
    phis = np.linspace(0.0, 2.0 * np.pi, M, endpoint=False)
    cos_phi = np.cos(phis)
    sin_phi = np.sin(phis)

    verts = []           # list of 3D vertex coordinates
    ring_indices = []    # for each profile index i, store list of M vertex indices (or repeated same index for degenerate ring)

    # Create vertices ring-by-ring. If radius <= tol_zero, create a single center vertex and
    # reference it M times to keep indexing consistent.
    for i in range(N):
        xi = float(x[i])
        ri = float(r[i])
        if ri <= tol_zero:
            # single center vertex
            idx_center = len(verts)
            verts.append((xi, 0.0, 0.0))
            ring_idx = [idx_center] * M
        else:
            ring_idx = []
            for j in range(M):
                y = ri * cos_phi[j]
                z = ri * sin_phi[j]
                verts.append((xi, y, z))
                ring_idx.append(len(verts) - 1)
        ring_indices.append(ring_idx)

    faces = []

    # Create side faces between adjacent rings (i -> i+1)
    for i in range(N - 1):
        A = ring_indices[i]
        B = ring_indices[i + 1]
        for j in range(M):
            j1 = (j + 1) % M
            a = A[j]
            b = A[j1]
            c = B[j1]
            d = B[j]
            # Add two triangles (a,b,c) and (a,c,d). Degenerate triangles are allowed
            # (they will be ignored/merged by trimesh.process if requested).
            faces.append((a, b, c))
            faces.append((a, c, d))

    # Caps: create fan triangles at start (i=0) and end (i=N-1) if requested and radius>tol_zero
    if cap_ends:
        # start cap (x = x[0]) -- orient so outside normal points in negative X (convention)
        if r[0] > tol_zero:
            # add center vertex
            idx_center_start = len(verts)
            verts.append((float(x[0]), 0.0, 0.0))
            A = ring_indices[0]
            for j in range(M):
                j1 = (j + 1) % M
                # triangle (center, A[j1], A[j]) gives normal pointing -X if rings ordered CCW looking down +X
                faces.append((idx_center_start, A[j1], A[j]))

        # end cap (x = x[-1]) -- orient so outside normal points in +X
        if r[-1] > tol_zero:
            idx_center_end = len(verts)
            verts.append((float(x[-1]), 0.0, 0.0))
            B = ring_indices[-1]
            for j in range(M):
                j1 = (j + 1) % M
                # triangle (center, B[j], B[j1]) gives normal pointing +X
                faces.append((idx_center_end, B[j], B[j1]))

    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    
    # =============================================================================
    print('np.shape(verts) =', np.shape(verts))
    print('np.shape(faces) =', np.shape(faces))
    verts[:,0] += 22
    # faces[:,0] += 30
    verts[:,1] += 10
    # faces[:,1] += 30
    print('np.shape(verts) =', np.shape(verts))
    print('np.shape(faces) =', np.shape(faces))
    # =============================================================================
    
    # print('np.shape(verts) =', np.shape(verts))
    # print('np.shape(faces) =', np.shape(faces))
    # sys.exit()

    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=process)

    # in the rare case processing fails to make watertight, attempt small fixes
    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        # as a last resort (if user accepts convex hull), return convex hull
        mesh = trimesh.Trimesh(vertices=verts, faces=faces).convex_hull

    return mesh

# -------------------------
# Demo & plotting
# -------------------------
if __name__ == '__main__':
    # Provided profile arrays (Psi_list_closed = x, zeta_list_closed = r)
    Psi_list_closed = np.array([
    0.   , 0.005, 0.01 , 0.015, 0.02 , 0.025, 0.03 , 0.035, 0.04 ,
    0.045, 0.05 , 0.055, 0.06 , 0.065, 0.07 , 0.075, 0.08 , 0.085,
    0.09 , 0.095, 0.1  , 0.105, 0.11 , 0.115, 0.12 , 0.125, 0.13 ,
    0.135, 0.14 , 0.145, 0.15 , 0.155, 0.16 , 0.165, 0.17 , 0.175,
    0.18 , 0.185, 0.19 , 0.195, 0.2  , 0.205, 0.21 , 0.215, 0.22 ,
    0.225, 0.23 , 0.235, 0.24 , 0.245, 0.25 , 0.255, 0.26 , 0.265,
    0.27 , 0.275, 0.28 , 0.285, 0.29 , 0.295, 0.3  , 0.305, 0.31 ,
    0.315, 0.32 , 0.325, 0.33 , 0.335, 0.34 , 0.345, 0.35 , 0.355,
    0.36 , 0.365, 0.37 , 0.375, 0.38 , 0.385, 0.39 , 0.395, 0.4  ,
    0.405, 0.41 , 0.415, 0.42 , 0.425, 0.43 , 0.435, 0.44 , 0.445,
    0.45 , 0.455, 0.46 , 0.465, 0.47 , 0.475, 0.48 , 0.485, 0.49 ,
    0.495, 0.5  , 0.505, 0.51 , 0.515, 0.52 , 0.525, 0.53 , 0.535,
    0.54 , 0.545, 0.55 , 0.555, 0.56 , 0.565, 0.57 , 0.575, 0.58 ,
    0.585, 0.59 , 0.595, 0.6  , 0.605, 0.61 , 0.615, 0.62 , 0.625,
    0.63 , 0.635, 0.64 , 0.645, 0.65 , 0.655, 0.66 , 0.665, 0.67 ,
    0.675, 0.68 , 0.685, 0.69 , 0.695, 0.7  , 0.705, 0.71 , 0.715,
    0.72 , 0.725, 0.73 , 0.735, 0.74 , 0.745, 0.75 , 0.755, 0.76 ,
    0.765, 0.77 , 0.775, 0.78 , 0.785, 0.79 , 0.795, 0.8  , 0.805,
    0.81 , 0.815, 0.82 , 0.825, 0.83 , 0.835, 0.84 , 0.845, 0.85 ,
    0.855, 0.86 , 0.865, 0.87 , 0.875, 0.88 , 0.885, 0.89 , 0.895,
    0.9  , 0.905, 0.91 , 0.915, 0.92 , 0.925, 0.93 , 0.935, 0.94 ,
    0.945, 0.95 , 0.955, 0.96 , 0.965, 0.97 , 0.975, 0.98 , 0.985,
    0.99 , 0.995, 1.   ]) * 4

    zeta_list_closed = np.array([
    0.        , 0.02741722, 0.03861335, 0.04705651, 0.05404363,
    0.06008184, 0.06543369, 0.07025658, 0.07465406, 0.07869873,
    0.08244396, 0.08593042, 0.08919002, 0.09224839, 0.09512649,
    0.09784178, 0.100409  , 0.1028407 , 0.10514776, 0.10733962,
    0.10942458, 0.11140999, 0.11330238, 0.1151076 , 0.11683093,
    0.11847712, 0.12005051, 0.12155505, 0.12299435, 0.12437173,
    0.12569026, 0.12695276, 0.12816185, 0.12931996, 0.13042937,
    0.13149217, 0.13251034, 0.13348573, 0.13442008, 0.13531502,
    0.13617206, 0.13699266, 0.13777818, 0.13852989, 0.13924902,
    0.13993671, 0.14059404, 0.14122206, 0.14182173, 0.14239398,
    0.1429397 , 0.14345972, 0.14395484, 0.14442581, 0.14487336,
    0.14529818, 0.14570091, 0.14608218, 0.14644258, 0.14678269,
    0.14710303, 0.14740412, 0.14768645, 0.14795049, 0.14819669,
    0.14842546, 0.14863722, 0.14883234, 0.14901119, 0.14917413,
    0.14932148, 0.14945357, 0.14957068, 0.14967311, 0.14976113,
    0.14983498, 0.14989493, 0.14994119, 0.14997399, 0.14999353,
    0.15      , 0.14999312, 0.14997261, 0.14993862, 0.14989133,
    0.14983087, 0.14975739, 0.14967102, 0.14957188, 0.14946008,
    0.14933573, 0.14919894, 0.14904979, 0.14888836, 0.14871473,
    0.14852898, 0.14833115, 0.14812132, 0.14789952, 0.14766581,
    0.14742021, 0.14716277, 0.1468935 , 0.14661242, 0.14631955,
    0.1460149 , 0.14569846, 0.14537024, 0.14503022, 0.1446784 ,
    0.14431476, 0.14393927, 0.14355191, 0.14315264, 0.14274142,
    0.14231823, 0.141883  , 0.14143569, 0.14097625, 0.14050461,
    0.14002072, 0.13952451, 0.1390159 , 0.13849482, 0.1379612 ,
    0.13741495, 0.13685598, 0.13628421, 0.13569953, 0.13510187,
    0.1344911 , 0.13386714, 0.13322987, 0.13257918, 0.13191495,
    0.13123708, 0.13054544, 0.1298399 , 0.12912034, 0.12838663,
    0.12763864, 0.12687623, 0.12609925, 0.12530758, 0.12450106,
    0.12367954, 0.12284288, 0.12199092, 0.12112351, 0.12024048,
    0.11934168, 0.11842693, 0.11749608, 0.11654896, 0.11558538,
    0.11460518, 0.11360818, 0.1125942 , 0.11156306, 0.11051456,
    0.10944854, 0.10836478, 0.10726311, 0.10614333, 0.10500524,
    0.10384864, 0.10267333, 0.10147911, 0.10026577, 0.0990331 ,
    0.0977809 , 0.09650894, 0.09521702, 0.09390492, 0.09257242,
    0.09121929, 0.08984532, 0.08845027, 0.08703393, 0.08559606,
    0.08413643, 0.0826548 , 0.08115095, 0.07962463, 0.0780756 ,
    0.07650362, 0.07490845, 0.07328985, 0.07164756, 0.06998133,
    0.06829092, 0.06657607, 0.06483652, 0.06307203, 0.06128233,
    0.05946715, 0.05762625, 0.05575935, 0.05386618, 0.05194649,
    0.05      ]) * 4

    # build mesh
    mesh = revolve_profile_to_mesh(Psi_list_closed, zeta_list_closed, n_phi=128//4, cap_ends=True)

    print("Mesh: verts", mesh.vertices.shape, "faces", mesh.faces.shape, "watertight?", mesh.is_watertight)

    # plot for verification
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection='3d')

    tri_verts = mesh.vertices[mesh.faces]
    mesh_collection = Poly3DCollection(tri_verts, alpha=0.85)
    mesh_collection.set_facecolor((0.6,0.8,1.0,0.9))
    mesh_collection.set_edgecolor((0.12,0.12,0.12,0.12))
    ax.add_collection3d(mesh_collection)

    # also plot profile curve (in the X-R plane)
    ax.plot(Psi_list_closed, zeta_list_closed, 0.0 * Psi_list_closed, 'k.-', label='profile (x,r)')

    # autoscale
    all_pts = mesh.vertices
    mins = all_pts.min(axis=0); maxs = all_pts.max(axis=0)
    ax.set_xlim(mins[0]-0.02, maxs[0]+0.02)
    ax.set_ylim(mins[1]-0.02, maxs[1]+0.02)
    ax.set_zlim(mins[2]-0.02, maxs[2]+0.02)
    ax.set_xlabel('X (Psi)'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title('Body of revolution from (Psi, zeta) profile')
    plt.tight_layout()
    ax.set_aspect('equal')
    plt.show()
