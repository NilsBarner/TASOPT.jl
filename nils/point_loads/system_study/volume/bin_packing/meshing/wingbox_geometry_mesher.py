__all__ = ["create_mesh_from_coordinates"]

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