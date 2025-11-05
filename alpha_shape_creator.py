__all__ = []

import sys
import alphashape
import numpy as np
import shapely.geometry as geom
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from shapely.geometry import Polygon, MultiPolygon


def add_shapely_polygon(ax, geom, facecolor='C0', edgecolor='k', alpha=0.3, zorder=1):
    """
    Safely add a shapely Polygon or MultiPolygon to a Matplotlib axis.
    Converts shapely coord sequences to explicit numpy arrays first to avoid
    errors from shapely/descartes versions.
    """
    if geom is None:
        return
    if getattr(geom, "is_empty", False):
        return

    def _add_ring(coords, **kwargs):
        # coords is an iterable of (x,y) tuples — convert to Nx2 array explicitly
        arr = np.asarray(list(coords))
        if arr.ndim != 2 or arr.shape[1] != 2 or arr.shape[0] < 3:
            return  # not a valid polygon ring
        patch = MplPolygon(arr, closed=True, **kwargs)
        ax.add_patch(patch)

    if isinstance(geom, Polygon):
        # exterior ring
        _add_ring(geom.exterior.coords, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, zorder=zorder)
        # interiors (holes): draw filled with the axis facecolor to cut them out visually
        bg = tuple(ax.get_facecolor())
        for interior in geom.interiors:
            _add_ring(interior.coords, facecolor=bg, edgecolor=None, alpha=1.0, zorder=zorder+1)
    elif isinstance(geom, MultiPolygon):
        for poly in geom.geoms:
            add_shapely_polygon(ax, poly, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, zorder=zorder)
    else:
        # Fallback: try to use exterior if present
        if hasattr(geom, "exterior"):
            _add_ring(geom.exterior.coords, facecolor=facecolor, edgecolor=edgecolor, alpha=alpha, zorder=zorder)


def ensure_xy_array(points):
    """
    Return an (N,2) numpy array of (x,y) from a variety of inputs:
    - list of (x,y) tuples
    - numpy array (N,2)
    - shapely.geometry.MultiPoint
    - shapely.geometry.Point
    - shapely.geometry.Polygon (use exterior coords)
    """
    if isinstance(points, np.ndarray):
        pts = points
    elif isinstance(points, list):
        pts = np.asarray(points)
    elif isinstance(points, geom.MultiPoint):
        pts = np.array([[p.x, p.y] for p in points.geoms])
    elif isinstance(points, geom.Point):
        pts = np.array([[points.x, points.y]])
    elif isinstance(points, geom.Polygon):
        pts = np.array(list(points.exterior.coords))
    else:
        # try a generic conversion
        try:
            pts = np.asarray(list(points))
        except Exception:
            raise TypeError(f"Can't convert points of type {type(points)} to (N,2) array")
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("Converted points must have shape (N,2)")
    return pts

#%%

if __name__ == '__main__':

    points_2d = np.load('points_2d.npy')
    pts = ensure_xy_array(points_2d)  # now pts is shape (N,2)
    
    # alpha_shape = alphashape.alphashape(pts, 0.00000005)
    alpha_shape = alphashape.alphashape(pts, 0.01)
    
    fig, ax = plt.subplots()
    ax.scatter(*zip(*points_2d), marker='.')
    add_shapely_polygon(
        ax, alpha_shape, facecolor='orange', edgecolor='None', alpha=0.3
    )
    ax.autoscale_view()
    plt.show()
