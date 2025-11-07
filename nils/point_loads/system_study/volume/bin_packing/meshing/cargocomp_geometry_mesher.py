#!/usr/bin/env python3
"""
cylindrical_double_chord_segment.py

Create a watertight mesh for the prism whose circular cross-section is the region
between two parallel chords of a circle. The chords are specified by angles
theta1 and theta2 measured from the positive X-axis (i.e. the circle diameter).
Both angles are taken in radians; we assume 0 <= theta1 <= theta2 <= pi/2 typically.

Inputs:
    r        : circle radius
    theta1   : angle (radians) between diameter and first chord (closer to +X if smaller)
    theta2   : angle (radians) between diameter and second chord (>= theta1)
    height   : extrusion distance along circle normal (positive by default)
    n_arc    : number of samples along each arc (>=2)

Returns:
    mesh, bottom_polygon_coords, top_polygon_coords

Visual verification plot is produced when run as __main__.
"""

__all__ = [
    "build_double_chord_polygon",
    "make_prism_from_polygon"
]

import sys
import numpy as np
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def build_double_chord_polygon(r, theta1, theta2, n_arc=64):
    """
    Build 2D coordinates (in XY plane, z=0) of polygon that represents the region
    between two parallel chords defined by angles theta1 and theta2.
    Ordering: top arc (theta1 -> theta2), chord at theta2 down, lower arc (-theta2 -> -theta1),
    chord at theta1 up.
    """
    if theta1 < 0 or theta2 < 0:
        raise ValueError("theta1, theta2 should be >= 0 (angles from +X diameter).")
    if theta2 < theta1:
        theta1, theta2 = theta2, theta1  # swap so theta1 <= theta2

    # sample top arc from theta1 to theta2 (inclusive)
    if n_arc < 2:
        n_arc = 2
    top_thetas = np.linspace(theta1, theta2, n_arc)
    top_arc = np.stack([-r * np.cos(top_thetas), r * np.sin(top_thetas), np.zeros_like(top_thetas)], axis=1)

    # chord at theta2: top point then bottom point
    x2 = -r * np.cos(theta2)
    y2 = r * np.sin(theta2)
    chord2 = np.array([[x2,  y2, 0.0],
                       [x2, -y2, 0.0]])

    # lower arc from -theta2 to -theta1 (inclusive)
    low_thetas = np.linspace(-theta2, -theta1, n_arc)
    low_arc = np.stack([-r * np.cos(low_thetas), r * np.sin(low_thetas), np.zeros_like(low_thetas)], axis=1)

    # chord at theta1: bottom then top (to close polygon)
    x1 = -r * np.cos(theta1)
    y1 = r * np.sin(theta1)
    chord1 = np.array([[x1, -y1, 0.0],
                       [x1,  y1, 0.0]])

    # assemble polygon: top_arc, chord2, low_arc, chord1
    poly = np.vstack([top_arc, chord2, low_arc, chord1])
    # print(np.shape(poly))
    # sys.exit()
    # =============================================================================
    poly[:,2] += 15
    # =============================================================================
    # remove consecutive duplicate points (if theta1==theta2 degenerate)
    # keep first occurrence
    uniq = [poly[0]]
    for p in poly[1:]:
        if np.linalg.norm(p - uniq[-1]) > 1e-12:
            uniq.append(p)
    poly = np.array(uniq)
    return poly


def make_prism_from_polygon(polygon_xy, height, extrude_direction=+1.0):
    """
    Given a 2D polygon (Nx3 with z=0), extrude it by 'height' along normal direction.
    extrude_direction = +1 or -1 chooses sign of extrusion (which side is "top").
    Returns watertight trimesh and bottom/top polygons.
    """
    bottom = polygon_xy.copy()
    top = bottom + np.array([0.0, 0.0, extrude_direction * height])

    Nb = bottom.shape[0]
    # build vertices: bottom [0..Nb-1], top [Nb..2Nb-1], then two centroids
    verts = np.vstack([bottom, top])
    centroid_bottom = bottom.mean(axis=0)
    centroid_top = top.mean(axis=0)
    idx_centroid_bottom = len(verts); verts = np.vstack([verts, centroid_bottom])
    idx_centroid_top = len(verts); verts = np.vstack([verts, centroid_top])

    faces = []

    # bottom cap triangles (ensure orientation pointing -extrude_direction)
    # Use fan triangulation about centroid_bottom.
    for i in range(Nb):
        i_next = (i + 1) % Nb
        # for bottom, we want normal to point opposite to extrude_direction.
        # use order (i_next, i, centroid) to get consistent orientation with typical right-hand rule
        faces.append([i_next, i, idx_centroid_bottom])

    # top cap triangles (orientation reversed)
    for i in range(Nb):
        i_next = (i + 1) % Nb
        ti = Nb + i; ti_next = Nb + i_next
        faces.append([ti, ti_next, idx_centroid_top])

    # sides: quads split into two triangles
    for i in range(Nb):
        i_next = (i + 1) % Nb
        b0 = i; b1 = i_next
        t0 = Nb + i; t1 = Nb + i_next
        faces.append([b0, b1, t1])
        faces.append([b0, t1, t0])

    verts = np.asarray(verts)
    faces = np.asarray(faces, dtype=np.int64)

    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=True)
    # try to fix minor non-watertightness
    if not mesh.is_watertight:
        try:
            mesh.remove_unreferenced_vertices()
            mesh.fill_holes()
        except Exception:
            pass
    if not mesh.is_watertight:
        mesh = trimesh.Trimesh(vertices=verts, faces=faces).convex_hull
        
    def swap_xz(mesh):
        new_vertices = mesh.vertices.copy()
        new_vertices[:, [0, 2]] = new_vertices[:, [2, 0]]
        mesh = trimesh.Trimesh(vertices=new_vertices, faces=mesh.faces, process=False)
        return mesh

    mesh = swap_xz(mesh)
    bottom[:, [0, 2]] = bottom[:, [2, 0]]
    top[:, [0, 2]] = top[:, [2, 0]]
        
    return mesh, bottom, top

# -------------------------
# Demo / verification plot
# -------------------------
if __name__ == "__main__":
    # Example parameters
    r = 3.0
    theta1_deg = 0.0
    thetafcs_deg = 55
    pidiv2minustheta1_deg = 90. - thetafcs_deg   # angle between diameter and first chord
    pidiv2minustheta2_deg = 90. - theta1_deg   # angle between diameter and second chord
    theta1 = np.deg2rad(pidiv2minustheta1_deg)
    theta2 = np.deg2rad(pidiv2minustheta2_deg)
    height = 15.2
    n_arc = 80

    # Build polygon and prism
    poly = build_double_chord_polygon(r, theta1, theta2, n_arc=n_arc)
    mesh, bottom, top = make_prism_from_polygon(poly, height, extrude_direction=+1.0)

    print("Polygon vertices (bottom) count:", bottom.shape[0])
    print("Mesh vertices, faces:", mesh.vertices.shape, mesh.faces.shape)
    print("Watertight?", mesh.is_watertight)

    # Plot
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111, projection='3d')

    tri_verts = mesh.vertices[mesh.faces]
    mesh_collection = Poly3DCollection(tri_verts)
    # mesh_collection.set_facecolor('red')
    # mesh_collection.set_edgecolor('red')
    # mesh_collection.set_alpha(0.1)
    mesh_collection.set(
        facecolor='red',
        edgecolor='none',
        alpha=0.1,
    )
    ax.add_collection3d(mesh_collection)

    # plot bottom polygon points
    # ax.scatter(bottom[:,0], bottom[:,1], bottom[:,2], color='red', s=40, label='bottom polygon verts')
    # plot chord endpoints explicitly
    x1, y1 = r*np.cos(theta1), r*np.sin(theta1)
    x2, y2 = r*np.cos(theta2), r*np.sin(theta2)
    # ax.scatter([x1,x1], [y1,-y1], [0,0], color='magenta', s=60, label='chord1 endpoints')
    # ax.scatter([x2,x2], [y2,-y2], [0,0], color='green', s=60, label='chord2 endpoints')

    # # plot polygon edges
    # def plot_poly_edges(ax, poly, color):
    #     for i in range(poly.shape[0]):
    #         j = (i+1) % poly.shape[0]
    #         seg = np.vstack([poly[i], poly[j]])
    #         ax.plot(seg[:,0], seg[:,1], seg[:,2], color=color, linewidth=1.5)
    # plot_poly_edges(ax, bottom, 'k')
    
    # =============================================================================
    # plot polygon edges
    def plot_poly_edges(ax, poly, color):
        for i in range(poly.shape[0]):
            j = (i+1) % poly.shape[0]
            seg = np.vstack([poly[i], poly[j]])
            ax.plot(seg[:,0], seg[:,1], seg[:,2], color=color, linewidth=1.5)
    
    # bottom edges (existing)
    plot_poly_edges(ax, bottom, 'k')
    # top edges (new) — top came from make_prism_from_polygon(...)
    plot_poly_edges(ax, top, 'k')
    
    # connect only the four corner points (chord endpoints)
    Nb = bottom.shape[0]
    # n_arc is the number of samples per arc used when building the polygon
    # (ensure n_arc from earlier scope is available; in your script it is).
    c2_top = int(n_arc)
    c2_bot = int(n_arc - 1)
    c1_bot = int(Nb - 2)
    c1_top = int(Nb - 1)
    
    corner_indices = [c2_top, c2_bot, c1_bot, c1_top]
    # corner_indices = [c2_top, c1_bot, c1_top]
    
    for i in corner_indices:
        seg = np.vstack([bottom[i], top[i]])
        ax.plot(seg[:,0], seg[:,1], seg[:,2], color='k', linewidth=1.5)
    # =============================================================================

    # autoscale & labels
    all_pts = mesh.vertices
    mins = all_pts.min(axis=0); maxs = all_pts.max(axis=0)
    ax.set_xlim(mins[0]-0.2, maxs[0]+0.2)
    ax.set_ylim(mins[1]-0.2, maxs[1]+0.2)
    ax.set_zlim(mins[2]-0.2, maxs[2]+0.2)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title(f'Cylindrical slab between chords: theta1={theta1_deg}°, theta2={thetafcs_deg}°, r={r}, h={height}')
    ax.legend()
    ax.set_aspect('equal')
    ax.view_init(elev=135, azim=-70, roll=-75)
    plt.show()
