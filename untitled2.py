import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.signal import savgol_filter  # for lightweight line smoothing

# --- example discrete 10x10 labels (replace with your data) ---
np.random.seed(0)
labels = np.random.randint(0, 4, (10, 10))

# --- colormap & norm for integer labels 0..3 ---
colors = ['#206095', '#a8bd3a', '#871a5b', '#f66068']
cmap = ListedColormap(colors)
norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5, 1.0), cmap.N)

fig, ax = plt.subplots(figsize=(6, 6))
# show raw cells (no interpolation of data)
ax.imshow(labels, origin='lower', interpolation='nearest', cmap=cmap, norm=norm,
          extent=[-0.5, labels.shape[1]-0.5, -0.5, labels.shape[0]-0.5])
ax.set_aspect('equal')

# compute contour lines at the half-integers (boundaries between labels)
levels = np.arange(0.5, len(colors) - 0.5 + 1)
CS = ax.contour(np.arange(labels.shape[1]), np.arange(labels.shape[0]),
                labels, levels=levels, linewidths=1.2, colors='k')

# robustly iterate over contour segments and smooth their vertex coordinates
for level_segs in CS.allsegs:            # each element is a list of segments for that level
    for seg in level_segs:                # seg is an (N,2) array of xy vertices
        if seg.shape[0] < 5:
            ax.plot(seg[:, 0], seg[:, 1], color='k', linewidth=1.2)
            continue
        # choose odd window length <= seg_length and >= 5 for savgol
        L = seg.shape[0]
        window = min(L if L % 2 == 1 else L-1, 51)  # cap window to 51 for speed
        if window < 5:
            ax.plot(seg[:, 0], seg[:, 1], color='k', linewidth=1.2)
            continue
        # smooth x and y independently (visual smoothing only)
        x_s = savgol_filter(seg[:, 0], window_length=window, polyorder=3, mode='interp')
        y_s = savgol_filter(seg[:, 1], window_length=window, polyorder=3, mode='interp')
        ax.plot(x_s, y_s, color='k', linewidth=1.2)

# legend / colorbar (optional)
from matplotlib.patches import Patch
legend_elems = [Patch(facecolor=colors[i], label=f'label {i}') for i in range(len(colors))]
ax.legend(handles=legend_elems, loc='upper right')

ax.set_xlim(-0.5, labels.shape[1]-0.5)
ax.set_ylim(-0.5, labels.shape[0]-0.5)
ax.set_xlabel('x'); ax.set_ylabel('y')
plt.tight_layout()
plt.show()
