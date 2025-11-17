import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.ndimage import zoom

# --- example discrete 10x10 data (labels 0..3) ---
np.random.seed(1)
labels = np.random.randint(0, 4, (10, 10))   # replace with your data

# --- upsample (smooth) the label field to a finer grid ---
upsample_factor = 10              # larger -> smoother boundaries
# zoom with spline interpolation (order=3) -> smooth continuous field
smooth_field = zoom(labels.astype(float), zoom=upsample_factor, order=3)

# --- coordinates for plotting ---
ny, nx = smooth_field.shape
x = np.linspace(0, labels.shape[1] - 1, nx)
y = np.linspace(0, labels.shape[0] - 1, ny)
X, Y = np.meshgrid(x, y)

# --- discrete colormap and boundaries for integer labels 0..3 ---
colors = ['#206095', '#a8bd3a', '#871a5b', '#f66068']
cmap = ListedColormap(colors)
n_levels = len(colors)
# boundaries positioned between integer label values so each color occupies [k-0.5, k+0.5)
bounds = np.arange(-0.5, n_levels + 0.5, 1.0)
norm = BoundaryNorm(bounds, cmap.N)

# --- plot: smooth-looking discrete regions via contourf ---
fig, ax = plt.subplots(figsize=(6, 5))
cf = ax.contourf(X, Y, smooth_field, levels=bounds, cmap=cmap, norm=norm, extend='neither')
ax.set_title("Smooth discrete contour plot (upsampled & interpolated)")
ax.set_xlabel("x")
ax.set_ylabel("y")

# colorbar with discrete ticks/labels
cbar = fig.colorbar(cf, ticks=np.arange(n_levels))
cbar.ax.set_yticklabels([f"label {i}" for i in range(n_levels)])

plt.tight_layout()
plt.show()
